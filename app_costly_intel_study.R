library(shiny)
library(shinyjs)
library(httr)
library(jsonlite)

# bayesian coverage allocator — study edition

# constants
N_SILOS        <- 10
START_UNITS    <- 9
TOTAL_MISSILES <- 8
RISK_K         <- 0.001   # very low: press6≈5%, press7≈8%, press8≈12% — comfortable 6-8 presses
RISK_EXP       <- 2.5     # steep acceleration keeps late presses meaningful
MAX_SAMPLES    <- 8       # max combined sample presses per trial
N_TRIALS       <- 4       # real trials (+ 1 practice = 5 total)
ADMIN_CODE     <- "admin"

# scripted hit/miss streams for fixed graph mode
# each list is one trial's sequence of hits (1) and misses (0)
SCRIPTED_STREAMS <- list(
  # Trial 0 (practice): A heavy, B light — obvious split
  list(A = c(1,1,0,1,1,0,1,1), B = c(0,0,1,0,0,1,0,0)),
  # Trial 1: A heavy — starts uncertain, converges high
  list(A = c(1,0,1,1,1,0,1,1), B = c(0,1,0,0,0,1,0,0)),
  # Trial 2: A heavy — quick convergence
  list(A = c(1,1,1,0,1,1,1,0), B = c(0,0,0,1,0,0,0,1)),
  # Trial 3: A heavy — slower convergence
  list(A = c(0,1,1,1,0,1,1,1), B = c(1,0,0,0,1,0,0,0)),
  # Trial 4: A heavy — very clear signal
  list(A = c(1,1,1,1,0,1,1,1), B = c(0,0,0,0,1,0,0,0))
)

# tutorial rounds shown before practice — probabilities revealed, no sampling
TUTORIAL_SCENARIOS <- list(
  list(
    label   = "Tutorial 1 of 3 — Easy",
    pA      = 0.80, pB = 0.20,
    mA      = 6,    mB = 2,
    hint    = "Field A is clearly more dangerous. Put most of your units there."
  ),
  list(
    label   = "Tutorial 2 of 3 — Moderate",
    pA      = 0.60, pB = 0.40,
    mA      = 5,    mB = 3,
    hint    = "Field A is somewhat more dangerous — but B still has missiles. Split carefully."
  ),
  list(
    label   = "Tutorial 3 of 3 — Tricky",
    pA      = 0.55, pB = 0.45,
    mA      = 5,    mB = 3,
    hint    = "The fields look similar. Even a small edge in probability matters here."
  )
)

# supabase connection
# fill these in before deploying
SUPABASE_URL      <- ""   # e.g. "https://abcdefgh.supabase.co"
SUPABASE_ANON_KEY <- ""   # your anon/public key (safe to expose in client)

# csv fallback if supabase is not configured
DATA_DIR <- "study_logs"
if (!dir.exists(DATA_DIR)) dir.create(DATA_DIR, recursive = TRUE)
EVENTS_FILE  <- file.path(DATA_DIR, "all_events.csv")
LAUNCHES_FILE <- file.path(DATA_DIR, "launch_trials.csv")

# participant codes — deterministic from seed 42

make_code_table <- function() {
  set.seed(42)
  n <- 50
  chars <- c(LETTERS, 0:9)
  codes <- replicate(n, paste0(sample(chars, 6, replace = TRUE), collapse = ""))
  # 10 codes per condition, 5 conditions
  modes <- rep(c("off", "curve", "dotplot", "gauge", "histogram"), each = 10)
  modes <- modes[sample(n)]   # shuffle (deterministically, seed already set)
  data.frame(code = codes, viz = modes, stringsAsFactors = FALSE)
}
CODE_TABLE <- make_code_table()
lookup_code <- function(code) {
  idx <- which(CODE_TABLE$code == toupper(code))
  if (length(idx) == 0) return(NULL)
  CODE_TABLE$viz[idx[1]]
}

# bayesian helpers
`%||%` <- function(a, b) if (!is.null(a)) a else b

detection_risk <- function(s, k, e) {
  if (s <= 0) return(0)
  1 - exp(-k * (s ^ e))
}

posterior_var <- function(a, b) (a * b) / (((a + b)^2) * (a + b + 1))

exp_covered <- function(a, b, n, k) {
  if (k <= 0) return(0)
  if (k >= n) return(n * a / (a + b))
  # Reduce grid for large n to keep computation fast
  grid_size <- min(500, max(100, 5000 %/% (n + 1)))
  p  <- seq(0.0001, 0.9999, length.out = grid_size)
  em <- sapply(p, function(pi) sum(pmin(0:n, k) * dbinom(0:n, n, pi)))
  sum(em * dbeta(p, a, b)) * (p[2] - p[1])
}

bayes_optimal_split <- function(a1, b1, a2, b2, n_silos, U) {
  if (U <= 0) return(list(opt_uA = 0, opt_uB = 0))
  # skip calculation if it would freeze the app
  if (U * n_silos > 5000) {
    return(list(opt_uA = NA, opt_uB = NA))
  }
  ka  <- 0:U
  eA  <- sapply(ka,     function(k) exp_covered(a1, b1, n_silos, k))
  eB  <- sapply(U - ka, function(k) exp_covered(a2, b2, n_silos, k))
  idx <- which.max(eA + eB)
  list(opt_uA = ka[idx], opt_uB = U - ka[idx])
}

make_scenarios <- function(n_total, total_missiles = TOTAL_MISSILES) {
  # field a is always heavier
  pA_pool <- c(0.70, 0.68, 0.72, 0.65, 0.70, 0.68, 0.72, 0.65)
  pB_pool <- c(0.30, 0.32, 0.28, 0.35, 0.30, 0.32, 0.28, 0.35)
  pA <- pA_pool[1:n_total]; pB <- pB_pool[1:n_total]
  set.seed(20260508)
  raw <- rbeta(n_total, 0.5, 0.5)
  # mA always gets strictly more than half; mB gets the rest
  # pmax/pmin ensures at least 1 missile each side
  mA  <- pmax(ceiling(total_missiles/2) + 1,
              pmin(total_missiles - 1, round(raw * (total_missiles - 2) + ceiling(total_missiles/2) + 1)))
  mB  <- total_missiles - mA
  lapply(seq_len(n_total), function(i) list(
    trial_index     = i - 1,
    pA_true         = pA[i], pB_true = pB[i],
    mA_true         = mA[i], mB_true = mB[i],
    sample_A_stream = rbinom(MAX_SAMPLES, 1, pA[i]),
    sample_B_stream = rbinom(MAX_SAMPLES, 1, pB[i])
  ))
}

# logging
supabase_insert <- function(table, row) {
  if (nchar(SUPABASE_URL) == 0 || nchar(SUPABASE_ANON_KEY) == 0) return(FALSE)
  tryCatch({
    resp <- httr::POST(
      paste0(SUPABASE_URL, "/rest/v1/", table),
      httr::add_headers(
        "apikey"        = SUPABASE_ANON_KEY,
        "Authorization" = paste("Bearer", SUPABASE_ANON_KEY),
        "Content-Type"  = "application/json",
        "Prefer"        = "return=minimal"
      ),
      body   = toJSON(list(row), auto_unbox = TRUE),
      encode = "raw"
    )
    httr::status_code(resp) %in% c(200, 201)
  }, error = function(e) FALSE)
}

append_csv <- function(row, file) {
  write.table(as.data.frame(row, stringsAsFactors = FALSE),
              file = file, sep = ",", row.names = FALSE,
              col.names = !file.exists(file),
              append = file.exists(file), qmethod = "double")
}

log_row <- function(row, table, csv_file) {
  if (!supabase_insert(table, row)) append_csv(row, csv_file)
}

# ui
ui <- fluidPage(
  useShinyjs(),
  tags$head(
    tags$link(rel = "preconnect", href = "https://fonts.googleapis.com"),
    tags$link(rel = "stylesheet",
      href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500;600&display=swap"),

    tags$script(HTML("
(function(){
var rnd=function(a,b){return a+Math.random()*(b-a);};
var W,H,AurC,AurX,ShpC,ShpX,AmbC,AmbX,SpkC,SpkX;

/* ── Aurora orbs ── */
var orbs=[],N_ORBS=6;
function mkOrb(){return{x:rnd(.05,.95),y:rnd(.05,.95),vx:rnd(-.00012,.00012),vy:rnd(-.00010,.00010),r:rnd(.20,.42),ph:rnd(0,6.28),sp:rnd(.004,.010),hue:rnd(135,172),sat:rnd(38,68),al:rnd(.030,.058)};}
function drawAurora(){AurX.clearRect(0,0,W,H);var mx=Math.max(W,H);for(var i=0;i<orbs.length;i++){var o=orbs[i];o.ph+=o.sp;o.x+=o.vx;o.y+=o.vy;if(o.x<-.2)o.x=1.1;if(o.x>1.2)o.x=-.1;if(o.y<-.2)o.y=1.1;if(o.y>1.2)o.y=-.1;var p=.55+.45*Math.sin(o.ph),cx=o.x*W,cy=o.y*H,rad=o.r*mx*p,al=o.al*p;var g=AurX.createRadialGradient(cx,cy,0,cx,cy,rad);g.addColorStop(0,'hsla('+o.hue+','+o.sat+'%,38%,'+al.toFixed(4)+')');g.addColorStop(.4,'hsla('+o.hue+','+o.sat+'%,25%,'+(al*.45).toFixed(4)+')');g.addColorStop(1,'hsla('+o.hue+','+o.sat+'%,12%,0)');AurX.beginPath();AurX.arc(cx,cy,rad,0,6.2832);AurX.fillStyle=g;AurX.fill();}requestAnimationFrame(drawAurora);}

/* ── Floating geometric shapes ── */
var shapes=[],NUM_SHAPES=14;
function mkShape(){var t=['tri','hex','sq','ring','pent'];return{x:rnd(0,W),y:rnd(0,H),vx:rnd(-.08,.08),vy:rnd(-.06,.06),rot:rnd(0,6.28),rs:rnd(-.003,.003),size:rnd(18,70),op:0,mop:rnd(.045,.13),fi:true,ir:rnd(.002,.005),or:rnd(.0005,.0015),type:t[Math.floor(Math.random()*t.length)],hue:rnd(130,195)};}
function polyPath(ctx,cx,cy,n,r,rot){ctx.beginPath();for(var i=0;i<=n;i++){var a=rot+i*(6.2832/n);i===0?ctx.moveTo(cx+r*Math.cos(a),cy+r*Math.sin(a)):ctx.lineTo(cx+r*Math.cos(a),cy+r*Math.sin(a));}ctx.closePath();}
function drawShapes(){ShpX.clearRect(0,0,W,H);for(var i=0;i<shapes.length;i++){var s=shapes[i];s.x+=s.vx;s.y+=s.vy;s.rot+=s.rs;s.fi?(s.op+=s.ir,s.op>=s.mop&&(s.op=s.mop,s.fi=false)):(s.op-=s.or);if(s.op<=0){shapes[i]=mkShape();continue;}if(s.x<-100)s.x=W+80;if(s.x>W+100)s.x=-80;if(s.y<-100)s.y=H+80;if(s.y>H+100)s.y=-80;ShpX.strokeStyle='hsla('+s.hue+',55%,62%,'+s.op.toFixed(3)+')';ShpX.lineWidth=0.9;if(s.type==='tri')polyPath(ShpX,s.x,s.y,3,s.size,s.rot);else if(s.type==='sq')polyPath(ShpX,s.x,s.y,4,s.size,s.rot);else if(s.type==='pent')polyPath(ShpX,s.x,s.y,5,s.size,s.rot);else if(s.type==='hex')polyPath(ShpX,s.x,s.y,6,s.size,s.rot);else{ShpX.beginPath();ShpX.arc(s.x,s.y,s.size,0,6.2832);ShpX.stroke();ShpX.beginPath();ShpX.arc(s.x,s.y,s.size*.62,0,6.2832);}ShpX.stroke();}requestAnimationFrame(drawShapes);}

/* ── Constellation ambient particles ── */
var pts=[],N=95,ACOLS=[[60,185,110],[45,162,90],[80,205,130],[55,172,100],[72,195,120],[48,158,85]];
function mkAmb(){var c=ACOLS[Math.floor(Math.random()*ACOLS.length)];return{x:rnd(0,W),y:rnd(0,H),vx:rnd(-.10,.10),vy:rnd(-.20,-.05),r:rnd(1.4,3.6),op:0,mop:rnd(.24,.48),fi:true,ir:rnd(.004,.009),or:rnd(.001,.003),r0:c[0],g0:c[1],b0:c[2]};}
function ambLoop(){var nw=window.innerWidth,nh=window.innerHeight;if(AmbC.width!==nw||AmbC.height!==nh){AmbC.width=nw;AmbC.height=nh;}W=nw;H=nh;AmbX.clearRect(0,0,W,H);
for(var i=0;i<pts.length;i++){var p=pts[i];p.x+=p.vx;p.y+=p.vy;p.fi?(p.op+=p.ir,p.op>=p.mop&&(p.op=p.mop,p.fi=false)):(p.op-=p.or);if(p.op<=0||p.y<-10){pts[i]=mkAmb();}}
var LINK=160;
for(var a=0;a<pts.length;a++){var pa=pts[a];if(pa.op<=0.06)continue;for(var b=a+1;b<pts.length;b++){var pb=pts[b];if(pb.op<=0.06)continue;var dx=pa.x-pb.x,dy=pa.y-pb.y,dist=Math.sqrt(dx*dx+dy*dy);if(dist<LINK){var la=(1-dist/LINK)*Math.min(pa.op,pb.op)*1.1;if(la<0.02)continue;AmbX.beginPath();AmbX.moveTo(pa.x,pa.y);AmbX.lineTo(pb.x,pb.y);AmbX.strokeStyle='rgba(90,210,140,'+Math.min(la,0.28).toFixed(3)+')';AmbX.lineWidth=1.0;AmbX.stroke();}}}
for(var j=0;j<pts.length;j++){var q=pts[j];if(q.op<=0)continue;AmbX.beginPath();AmbX.arc(q.x,q.y,q.r,0,6.2832);AmbX.fillStyle='rgba('+q.r0+','+q.g0+','+q.b0+','+q.op.toFixed(3)+')';AmbX.fill();}
requestAnimationFrame(ambLoop);}

/* ── Sparks + rings on click ── */
var sparks=[],rings=[],SCOLS=[[85,210,135],[65,175,105],[105,225,155],[125,235,165],[72,195,115]];
function mkSpark(x,y){var c=SCOLS[Math.floor(Math.random()*SCOLS.length)];var ang=Math.random()*6.2832,sp=rnd(1.2,5.0);return{x:x,y:y,vx:Math.cos(ang)*sp,vy:Math.sin(ang)*sp,r:rnd(1.2,3.2),op:rnd(.65,.92),dc:rnd(.022,.050),r0:c[0],g0:c[1],b0:c[2]};}
function mkRing(x,y){return{x:x,y:y,r:2,mr:rnd(40,90),op:.55,dc:.022};}
function sparkLoop(){SpkX.clearRect(0,0,W,H);
for(var j=rings.length-1;j>=0;j--){var rg=rings[j];rg.r+=(rg.mr-rg.r)*.09;rg.op-=rg.dc;if(rg.op<=0){rings.splice(j,1);continue;}SpkX.beginPath();SpkX.arc(rg.x,rg.y,rg.r,0,6.2832);SpkX.strokeStyle='rgba(100,220,150,'+rg.op.toFixed(3)+')';SpkX.lineWidth=1.5;SpkX.stroke();SpkX.beginPath();SpkX.arc(rg.x,rg.y,rg.r*.62,0,6.2832);SpkX.strokeStyle='rgba(60,190,110,'+(rg.op*.4).toFixed(3)+')';SpkX.lineWidth=.8;SpkX.stroke();}
for(var i=sparks.length-1;i>=0;i--){var p=sparks[i];p.x+=p.vx;p.y+=p.vy;p.vy+=.09;p.op-=p.dc;if(p.op<=0){sparks.splice(i,1);continue;}SpkX.beginPath();SpkX.arc(p.x,p.y,p.r,0,6.2832);SpkX.fillStyle='rgba('+p.r0+','+p.g0+','+p.b0+','+p.op.toFixed(3)+')';SpkX.fill();}
requestAnimationFrame(sparkLoop);}

/* ── Audio ── */
var AC=null;
function getAC(){if(!AC)try{AC=new(window.AudioContext||window.webkitAudioContext)();}catch(e){}return AC;}
function playClick(){var ac=getAC();if(!ac)return;var o=ac.createOscillator(),g=ac.createGain();o.connect(g);g.connect(ac.destination);o.type='sine';o.frequency.setValueAtTime(420,ac.currentTime);o.frequency.exponentialRampToValueAtTime(260,ac.currentTime+.06);g.gain.setValueAtTime(.06,ac.currentTime);g.gain.exponentialRampToValueAtTime(.0001,ac.currentTime+.09);o.start(ac.currentTime);o.stop(ac.currentTime+.1);}
function playHover(){var ac=getAC();if(!ac)return;var o=ac.createOscillator(),g=ac.createGain();o.connect(g);g.connect(ac.destination);o.type='sine';o.frequency.setValueAtTime(520,ac.currentTime);o.frequency.exponentialRampToValueAtTime(480,ac.currentTime+.04);g.gain.setValueAtTime(.025,ac.currentTime);g.gain.exponentialRampToValueAtTime(.0001,ac.currentTime+.05);o.start(ac.currentTime);o.stop(ac.currentTime+.06);}
function attachSounds(){
  document.addEventListener('click',function(e){if(e.target.closest('button:not(:disabled)'))playClick();});
  document.addEventListener('mouseover',function(e){var btn=e.target.closest('button:not(:disabled)');if(btn&&btn!==btn._lh){btn._lh=btn;playHover();setTimeout(function(){btn._lh=null;},200);}});
}

/* ── Canvas setup + init ── */
function mkC(zi){var c=document.createElement('canvas');c.style.cssText='position:fixed;top:0;left:0;width:100%;height:100%;pointer-events:none;z-index:'+zi+';';c.width=W;c.height=H;document.body.appendChild(c);return c;}
function init(){W=window.innerWidth;H=window.innerHeight;
AurC=mkC(1);AurX=AurC.getContext('2d');
ShpC=mkC(1);ShpX=ShpC.getContext('2d');
AmbC=mkC(1);AmbX=AmbC.getContext('2d');
SpkC=mkC(999);SpkX=SpkC.getContext('2d');
for(var i=0;i<N_ORBS;i++)orbs.push(mkOrb());
for(var j=0;j<NUM_SHAPES;j++)shapes.push(mkShape());
for(var k=0;k<N;k++)pts.push(mkAmb());
drawAurora();drawShapes();ambLoop();sparkLoop();
document.addEventListener('click',function(e){var n=12+Math.floor(Math.random()*8);for(var i=0;i<n;i++)sparks.push(mkSpark(e.clientX,e.clientY));rings.push(mkRing(e.clientX,e.clientY));});
attachSounds();
window.addEventListener('resize',function(){W=window.innerWidth;H=window.innerHeight;[AurC,ShpC,SpkC].forEach(function(c){c.width=W;c.height=H;});});
}
document.readyState==='loading'?document.addEventListener('DOMContentLoaded',init):setTimeout(init,60);
(function(){
var pA1=1,pB1=1,pA2=1,pB2=1,cA1=1,cB1=1,cA2=1,cB2=1,raf=null,ov=null,OX=null;
function bpdf(x,a,b){return(a<=0||b<=0)?1:Math.pow(x,a-1)*Math.pow(1-x,b-1);}
function peak(a,b){return(a>1&&b>1)?bpdf((a-1)/(a+b-2),a,b):(bpdf(.5,a,b)||1);}
function draw(t){if(!ov||!OX)return;var W=ov.width,H=ov.height;OX.clearRect(0,0,W,H);var e=1-Math.pow(1-t,3);var a1=pA1+(cA1-pA1)*e,b1=pB1+(cB1-pB1)*e,a2=pA2+(cA2-pA2)*e,b2=pB2+(cB2-pB2)*e;var pkA=peak(a1,b1)||1,pkB=peak(a2,b2)||1,ys=H*.76/Math.max(pkA,pkB,.5),steps=220;OX.beginPath();for(var i=0;i<=steps;i++){var x=.001+i/steps*.998,y=H-(bpdf(x,a1,b1)/Math.max(pkA,.01))*ys*.88;i?OX.lineTo(x*W,y):OX.moveTo(x*W,y);}OX.strokeStyle='rgba(201,96,112,'+(.45+.55*t).toFixed(2)+')';OX.lineWidth=2.5;OX.stroke();OX.beginPath();for(var j=0;j<=steps;j++){var xj=.001+j/steps*.998,yj=H-(bpdf(xj,a2,b2)/Math.max(pkB,.01))*ys*.88;j?OX.lineTo(xj*W,yj):OX.moveTo(xj*W,yj);}OX.strokeStyle='rgba(61,134,168,'+(.45+.55*t).toFixed(2)+')';OX.lineWidth=2.5;OX.stroke();ov.style.opacity=(1-t).toFixed(3);}
function tween(){if(raf)cancelAnimationFrame(raf);var dur=480,s=null;function step(ts){if(!s)s=ts;var p=Math.min((ts-s)/dur,1);draw(p);p<1?raf=requestAnimationFrame(step):ov&&(ov.style.opacity=0);}raf=requestAnimationFrame(step);}
function attachOv(){var cont=document.getElementById('plot');if(!cont){setTimeout(attachOv,500);return;}var par=cont.parentElement;ov=document.createElement('canvas');ov.style.cssText='position:absolute;top:0;left:0;width:100%;height:100%;pointer-events:none;z-index:3;opacity:0;';par.style.position='relative';par.appendChild(ov);function size(){var img=cont.querySelector('img');if(img&&img.offsetWidth){ov.width=img.offsetWidth;ov.height=img.offsetHeight;OX=ov.getContext('2d');}}var obs=new MutationObserver(function(){size();var d=document.getElementById('plotData');if(d){var na1=+d.dataset.a1||1,nb1=+d.dataset.b1||1,na2=+d.dataset.a2||1,nb2=+d.dataset.b2||1;if(na1!==cA1||nb1!==cB1||na2!==cA2||nb2!==cB2){pA1=cA1;pB1=cB1;pA2=cA2;pB2=cB2;cA1=na1;cB1=nb1;cA2=na2;cB2=nb2;tween();}}});obs.observe(cont,{childList:true,subtree:true,attributes:true});}
document.readyState==='loading'?document.addEventListener('DOMContentLoaded',attachOv):setTimeout(attachOv,600);
})();
})();
    ")),

    tags$style(HTML("
*,*::before,*::after{box-sizing:border-box;margin:0;padding:0;}
:root{
  --bg:#060f0a;--srf:#0a1610;--srf2:#0f1e16;--srf3:#162a1e;
  --bdr:#1a3024;--bdr2:#22402e;
  --tx:#d8efe2;--tx2:#7ec49a;--tx3:#4a8c68;
  --acc:#3d9962;--acc2:#2d7a4a;--accb:#061008;--accm:#143d22;
  --red:#c96070;--redb:#130509;--redm:#381018;
  --blue:#3d86a8;--blueb:#040d12;--bluem:#0d2230;
  --grn:#3db87a;--grnb:#040e08;--grnm:#0c3320;
  --amb:#b88a38;--ambb:#120900;--ambm:#382600;
  --sans:'Inter',-apple-system,sans-serif;
  --mono:'JetBrains Mono','Courier New',monospace;
  --r1:5px;--r2:9px;--r3:14px;
  --sh:0 2px 12px rgba(0,0,0,.55),0 0 0 1px rgba(255,255,255,.025);
  --sh2:0 6px 32px rgba(0,0,0,.70);
  --sh3:0 0 0 1px rgba(61,153,98,.25),0 4px 20px rgba(61,153,98,.08);
}
html{background:
  radial-gradient(ellipse 70% 45% at 15% 35%, rgba(45,122,74,.09) 0%, transparent 55%),
  radial-gradient(ellipse 55% 38% at 82% 68%, rgba(30,80,120,.07) 0%, transparent 50%),
  radial-gradient(ellipse 60% 50% at 52% 8%,  rgba(61,153,98,.05) 0%, transparent 45%),
  #060f0a !important;}
body{background:transparent!important;color:var(--tx)!important;font-family:var(--sans)!important;font-size:14px;line-height:1.55;-webkit-font-smoothing:antialiased;min-height:100vh;}
.container-fluid{position:relative;z-index:2;}
.modal-backdrop{z-index:1040!important;}.modal{z-index:1050!important;}

@keyframes fadeUp{from{opacity:0;transform:translateY(10px)}to{opacity:1;transform:none}}
@keyframes fadeIn{from{opacity:0}to{opacity:1}}
@keyframes bannerIn{from{opacity:0;transform:translateY(-6px) scale(.98)}to{opacity:1;transform:none}}
@keyframes shake{0%,100%{transform:translateX(0)}18%{transform:translateX(-5px)}36%{transform:translateX(5px)}54%{transform:translateX(-3px)}72%{transform:translateX(3px)}}
@keyframes pillPulse{0%,100%{box-shadow:0 0 0 0 rgba(61,153,98,0)}50%{box-shadow:0 0 0 5px rgba(61,153,98,.16)}}
@keyframes txtIn{from{opacity:0;transform:translateY(7px)}to{opacity:1;transform:none}}
@keyframes typeReveal{from{clip-path:inset(0 100% 0 0);opacity:.15}to{clip-path:inset(0 0% 0 0);opacity:1}}
@keyframes screenIn{from{opacity:0;transform:translateY(30px) scale(.97)}to{opacity:1;transform:none}}
@keyframes numPop{from{transform:translateY(4px) scale(.92);opacity:.4}to{transform:none;opacity:1}}
@keyframes shimmer{0%{background-position:-200% center}100%{background-position:200% center}}
@keyframes plotFade{from{opacity:.1;transform:scale(.994) translateY(5px)}to{opacity:1;transform:none}}
@keyframes vizIn{from{opacity:0;transform:translateY(8px)}to{opacity:1;transform:none}}
@keyframes cardIn{from{opacity:0;transform:translateY(12px) scale(.98)}to{opacity:1;transform:none}}
@keyframes rsStatIn{from{opacity:0;transform:translateY(16px)}to{opacity:1;transform:none}}
@keyframes rsRowIn{from{opacity:0;transform:translateX(-14px)}to{opacity:1;transform:none}}
@keyframes rsNumIn{from{opacity:0;transform:scale(.70)}to{opacity:1;transform:none}}
@keyframes glowPulse{0%,100%{box-shadow:var(--sh)}50%{box-shadow:var(--sh),0 0 22px rgba(61,153,98,.14)}}

.wrap{max-width:1280px;margin:0 auto;padding:18px 22px;display:flex;flex-direction:column;gap:12px;}
.banner-slot{min-height:66px;}

/* ── Topbar with shimmer ── */
.topbar{display:flex;justify-content:space-between;align-items:center;padding:12px 18px;
  background:var(--srf);border:1px solid var(--bdr2);border-top:2px solid var(--acc2);
  border-radius:var(--r3);box-shadow:var(--sh);
  background-image:linear-gradient(90deg,transparent 0%,rgba(61,153,98,.07) 40%,rgba(100,200,140,.12) 50%,rgba(61,153,98,.07) 60%,transparent 100%);
  background-size:200% 100%;
  animation:fadeUp .4s cubic-bezier(.22,1,.36,1),shimmer 8s linear 1.5s infinite!important;}
.topbar-brand{display:flex;align-items:center;gap:10px;}
.brand-icon{width:32px;height:32px;background:var(--accb);border:1px solid var(--accm);border-radius:8px;display:flex;align-items:center;justify-content:center;font-size:16px;transition:box-shadow .3s;}
.brand-icon:hover{box-shadow:0 0 12px rgba(61,153,98,.25);}
.topbar-title{font-size:11px;font-weight:700;letter-spacing:.18em;text-transform:uppercase;color:var(--tx3);}
.topbar-right{display:flex;align-items:center;gap:8px;}

/* ── Grid / Cards ── */
.main-grid{display:grid;grid-template-columns:358px 1fr;gap:12px;align-items:start;}
.left-col{display:flex;flex-direction:column;gap:12px;}
.card{background:var(--srf);border:1px solid var(--bdr);border-radius:var(--r3);padding:16px 18px;
  box-shadow:var(--sh);animation:cardIn .4s cubic-bezier(.22,1,.36,1) both;
  transition:box-shadow .25s ease,border-color .25s ease,transform .2s ease;}
.card:hover{box-shadow:var(--sh3),var(--sh)!important;border-color:rgba(61,153,98,.3)!important;transform:translateY(-1px);}
.left-col .card:nth-child(1){animation-delay:.05s}
.left-col .card:nth-child(2){animation-delay:.12s}
.left-col .card:nth-child(3){animation-delay:.19s}
.chart-card{background:var(--srf);border:1px solid var(--bdr);border-radius:var(--r3);padding:16px 18px;
  box-shadow:var(--sh);animation:cardIn .4s cubic-bezier(.22,1,.36,1) .08s both;
  min-height:500px;display:flex;flex-direction:column;
  transition:box-shadow .25s ease,border-color .25s ease,transform .2s ease;}
.chart-card:hover{box-shadow:var(--sh3),var(--sh)!important;border-color:rgba(61,153,98,.3)!important;transform:translateY(-1px);}
.card-title{font-size:10px;font-weight:700;letter-spacing:.14em;text-transform:uppercase;color:var(--tx3);margin-bottom:14px;}
.consent-screen{display:flex;align-items:center;justify-content:center;min-height:90vh;animation:screenIn .6s cubic-bezier(.22,1,.36,1) both;}
.consent-card{background:var(--srf);border:1px solid var(--bdr2);border-top:3px solid var(--acc);border-radius:var(--r3);padding:44px 48px 40px;max-width:640px;width:100%;box-shadow:var(--sh2);}
.consent-title{font-size:24px;font-weight:800;color:var(--tx);text-align:center;margin-bottom:8px;animation:typeReveal .7s cubic-bezier(.22,1,.36,1) .1s both;}
.consent-sub{font-family:var(--mono);font-size:11px;color:var(--tx3);text-align:center;letter-spacing:.1em;text-transform:uppercase;margin-bottom:28px;animation:txtIn .5s ease .3s both;}
.consent-divider{height:1px;background:var(--bdr);margin:22px 0;}
.consent-section{font-size:10px;font-weight:700;letter-spacing:.14em;text-transform:uppercase;color:var(--tx3);margin-bottom:10px;}
.consent-body{font-size:13px;color:#9ecfb5;line-height:1.85;margin-bottom:20px;}
.consent-body strong{color:var(--tx);}
.consent-code-row{display:flex;flex-direction:column;gap:6px;background:var(--srf2);border:1px solid var(--bdr);border-radius:var(--r2);padding:14px 16px;margin-bottom:14px;}
.consent-hint{font-size:11px;color:#9ecfb5;margin-top:3px;}
.consent-btns{display:flex;gap:10px;margin-top:4px;}
.consent-btns button{flex:1;}

/* ── Topbar ── */
.topbar{display:flex;justify-content:space-between;align-items:center;padding:12px 18px;background:var(--srf);border:1px solid var(--bdr2);border-top:2px solid var(--acc2);border-radius:var(--r3);box-shadow:var(--sh);animation:fadeUp .3s ease;}
.topbar-brand{display:flex;align-items:center;gap:10px;}
.brand-icon{width:32px;height:32px;background:var(--accb);border:1px solid var(--accm);border-radius:8px;display:flex;align-items:center;justify-content:center;font-size:16px;}
.topbar-title{font-size:11px;font-weight:700;letter-spacing:.18em;text-transform:uppercase;color:var(--tx3);}
.topbar-right{display:flex;align-items:center;gap:8px;}

/* ── Grid / Cards ── */
.main-grid{display:grid;grid-template-columns:358px 1fr;gap:12px;align-items:start;}
.left-col{display:flex;flex-direction:column;gap:12px;}
.card{background:var(--srf);border:1px solid var(--bdr);border-radius:var(--r3);padding:16px 18px;box-shadow:var(--sh);animation:fadeUp .35s ease both;}
.card:hover{box-shadow:var(--sh),0 0 0 1px rgba(61,153,98,.18),0 0 18px rgba(61,153,98,.08)!important;border-color:rgba(61,153,98,.35)!important;transition:all .25s!important;}
.chart-card{background:var(--srf);border:1px solid var(--bdr);border-radius:var(--r3);padding:16px 18px;box-shadow:var(--sh);animation:fadeUp .35s ease .06s both;min-height:500px;display:flex;flex-direction:column;}
.chart-card:hover{box-shadow:var(--sh),0 0 0 1px rgba(61,153,98,.18),0 0 18px rgba(61,153,98,.08)!important;border-color:rgba(61,153,98,.35)!important;transition:all .25s!important;}
.left-col .card:nth-child(1){animation-delay:.03s}.left-col .card:nth-child(2){animation-delay:.08s}.left-col .card:nth-child(3){animation-delay:.13s}
.card-title{font-size:10px;font-weight:700;letter-spacing:.14em;text-transform:uppercase;color:var(--tx3);margin-bottom:14px;}

/* ── Pill ── */
.pill{display:inline-flex;align-items:center;gap:5px;padding:5px 12px;border-radius:999px;font-family:var(--mono);font-size:11px;font-weight:500;border:1px solid var(--bdr2);background:var(--srf2);color:var(--tx2);white-space:nowrap;}
.pill.active{border-color:var(--acc);background:var(--accb);color:var(--acc);animation:pillPulse 2.5s ease-in-out infinite;}
.pill.practice-pill{border-color:var(--ambm);background:var(--ambb);color:var(--amb);}

/* ── Buttons ── */
button,.btn{font-family:var(--sans)!important;font-size:13px!important;font-weight:600!important;padding:8px 14px!important;border-radius:var(--r2)!important;border:1px solid var(--bdr2)!important;background:var(--srf2)!important;color:var(--tx)!important;cursor:pointer;line-height:1.4!important;transition:background .14s,border-color .14s,box-shadow .14s,transform .12s ease!important;}
button:hover:not(:disabled){background:var(--srf3)!important;border-color:var(--tx3)!important;box-shadow:0 4px 12px rgba(0,0,0,.45)!important;transform:translateY(-2px) scale(1.025)!important;}
button:active:not(:disabled){transform:translateY(0) scale(.98)!important;transition-duration:.06s!important;}
button:disabled{opacity:.3!important;cursor:default!important;}
@keyframes samplePulse{0%,100%{box-shadow:0 0 0 0 rgba(61,153,98,0),var(--sh)}60%{box-shadow:0 0 0 7px rgba(61,153,98,.13),var(--sh)}}
@keyframes launchGlow{0%,100%{box-shadow:var(--sh)}50%{box-shadow:0 0 16px rgba(61,153,98,.28),var(--sh)}}
@keyframes nextAmber{0%,100%{box-shadow:0 0 0 0 rgba(184,138,56,0)}60%{box-shadow:0 0 0 5px rgba(184,138,56,.14)}}
.btn-primary{background:var(--acc)!important;border-color:var(--acc)!important;color:#030a05!important;font-weight:700!important;animation:launchGlow 3s ease-in-out infinite;}
.btn-primary:hover:not(:disabled){background:#4aac72!important;border-color:#4aac72!important;box-shadow:0 0 20px rgba(61,153,98,.30)!important;animation:none!important;}
.btn-next{background:var(--ambb)!important;border-color:var(--ambm)!important;color:var(--amb)!important;font-weight:600!important;animation:nextAmber 3.4s ease-in-out infinite;}
.btn-next:hover:not(:disabled){background:var(--ambm)!important;animation:none!important;}
.btn-decline{background:var(--redb)!important;border-color:var(--redm)!important;color:var(--red)!important;font-weight:600!important;}
.btn-decline:hover:not(:disabled){background:var(--redm)!important;}
.btn-cancel{background:var(--srf3)!important;border-color:var(--bdr2)!important;color:var(--tx2)!important;}
.btn-cancel:hover:not(:disabled){background:#1c3428!important;color:var(--tx)!important;}
.topbar-btn{display:inline-flex!important;align-items:center!important;gap:5px!important;padding:6px 13px!important;border-radius:999px!important;font-family:var(--sans)!important;font-size:12px!important;font-weight:600!important;border:1px solid var(--bdr2)!important;background:transparent!important;color:var(--tx2)!important;transition:all .14s!important;}
.topbar-btn:hover:not(:disabled){background:var(--srf2)!important;border-color:var(--acc2)!important;color:var(--tx)!important;}
.topbar-btn.on{background:var(--accb)!important;border-color:var(--acc)!important;color:var(--acc)!important;}
.btn-sample{background:var(--srf2)!important;border:1px solid var(--acc2)!important;color:var(--tx)!important;font-size:14px!important;font-weight:700!important;padding:13px 18px!important;width:100%;margin-bottom:0;letter-spacing:.03em;animation:samplePulse 2.6s ease-in-out infinite;}
.btn-sample:hover:not(:disabled){background:var(--accb)!important;border-color:var(--acc)!important;color:var(--acc)!important;box-shadow:0 0 18px rgba(61,153,98,.20)!important;animation:none!important;}
.btn-sample:active:not(:disabled){transform:translateY(0) scale(.98)!important;}

/* ── Info strips ── */
.settings-display{font-family:var(--mono);font-size:12px;color:var(--tx2);background:var(--srf2);border:1px solid var(--bdr);border-radius:var(--r2);padding:10px 13px;margin-bottom:10px;line-height:1.9;}
.settings-display strong{color:var(--tx)}.badge-risk{color:#c87a4d;font-weight:700;}
.status-strip{display:flex;flex-wrap:wrap;gap:5px 14px;font-family:var(--mono);font-size:12px;background:var(--srf2);border:1px solid var(--bdr);border-radius:var(--r2);padding:8px 13px;}
.kv{display:flex;gap:4px;align-items:center}.k{color:var(--tx3)}.v{color:var(--tx);font-weight:600}.sep{color:var(--bdr2)}

/* ── Field rows ── */
.field-row{display:flex;align-items:center;gap:10px;border-radius:var(--r2);padding:10px 13px;margin-bottom:8px;border:1px solid transparent;}
.field-row.red{background:var(--redb);border-color:var(--redm)}.field-row.blue{background:var(--blueb);border-color:var(--bluem)}
.fr-tag{font-size:10px;font-weight:700;letter-spacing:.12em;text-transform:uppercase;min-width:52px}
.fr-tag.red{color:var(--red)}.fr-tag.blue{color:var(--blue)}
.fr-val{font-family:var(--mono);font-size:20px;font-weight:700;min-width:64px;line-height:1}
.fr-val.red{color:var(--red)}.fr-val.blue{color:var(--blue)}
.fr-meta{flex:1;font-family:var(--mono);font-size:11px;color:var(--tx2);display:flex;flex-direction:column;gap:3px}

/* ── Risk gauge ── */
.risk-gauge{background:var(--ambb);border:1px solid var(--ambm);border-radius:var(--r2);padding:11px 13px;margin-bottom:10px;transition:border-color .3s}
.risk-gauge.hi{border-color:#6a2020}
.risk-top{display:flex;justify-content:space-between;align-items:center;margin-bottom:8px}
.risk-title{font-size:10px;font-weight:700;letter-spacing:.14em;text-transform:uppercase;color:var(--amb)}
.risk-pct{font-family:var(--mono);font-size:24px;font-weight:700;line-height:1;transition:color .3s}
.risk-track{height:5px;background:rgba(184,138,56,.12);border-radius:3px;overflow:hidden;margin-bottom:7px}
.risk-fill{height:100%;border-radius:3px;transition:width .4s ease,background .4s ease}
.risk-hint{font-family:var(--mono);font-size:10px;color:var(--tx3);line-height:1.7}

/* ── Alloc slider ── */
.alloc-header{display:flex;justify-content:space-between;align-items:baseline;margin-bottom:6px}
.alloc-label{font-size:10px;font-weight:700;letter-spacing:.12em;text-transform:uppercase;color:var(--tx3)}
.alloc-b{font-family:var(--mono);font-size:18px;font-weight:700;color:var(--blue)}
.alloc-row{margin-bottom:12px}.alloc-row .control-label{display:none!important}
.irs--shiny .irs-line{background:var(--srf3);border-color:var(--bdr2);height:5px}
.irs--shiny .irs-bar{background:var(--acc);border:none;height:5px}
.irs--shiny .irs-handle{background:var(--acc);border:2px solid var(--srf);border-radius:50%;width:17px;height:17px;top:18px;box-shadow:0 0 5px rgba(61,153,98,.3)}
.irs--shiny .irs-single{background:var(--srf3);color:var(--tx2);font-family:var(--mono);font-size:11px;border:1px solid var(--bdr2);border-radius:4px}
/* Field A slider — coral */
#uA ~ .irs .irs-bar{background:#ff8096!important}
#uA ~ .irs .irs-handle{background:#ff8096!important;box-shadow:0 0 5px rgba(255,128,150,.3)!important}
/* Field B slider — blue */
#uB_slider ~ .irs .irs-bar{background:#5bc4f0!important}
#uB_slider ~ .irs .irs-handle{background:#5bc4f0!important;box-shadow:0 0 5px rgba(91,196,240,.3)!important}
/* Tutorial sliders */
#tut_uA ~ .irs .irs-bar{background:#ff8096!important}
#tut_uA ~ .irs .irs-handle{background:#ff8096!important;box-shadow:0 0 5px rgba(255,128,150,.3)!important}
#tut_uB ~ .irs .irs-bar{background:#5bc4f0!important}
#tut_uB ~ .irs .irs-handle{background:#5bc4f0!important;box-shadow:0 0 5px rgba(91,196,240,.3)!important}
.divider{height:1px;background:var(--bdr);margin:12px 0}

/* ── Banners ── */
.banner{border-radius:var(--r2);padding:12px 16px;display:flex;flex-direction:column;gap:5px;margin-bottom:4px}
.banner.win{background:var(--grnb);border:1px solid #1a5c34;animation:bannerIn .28s ease}
.banner.lose{background:var(--redb);border:1px solid #5a1a28;animation:bannerIn .28s ease,shake .38s ease .1s}
.banner.det{background:var(--ambb);border:1px solid #5a3a10;animation:bannerIn .28s ease,shake .38s ease .1s}
.banner.prac{background:var(--srf2);border:1px solid var(--bdr2);animation:bannerIn .28s ease}
.banner-tag{font-size:10px;font-weight:700;letter-spacing:.14em;text-transform:uppercase;display:block;animation:txtIn .32s ease .04s both}
.banner.win .banner-tag{color:var(--grn)}.banner.lose .banner-tag{color:var(--red)}.banner.det .banner-tag{color:var(--amb)}.banner.prac .banner-tag{color:var(--tx3)}
.banner-main{font-size:14px;font-weight:600;color:var(--tx);display:block;animation:txtIn .34s ease .10s both}
.banner-sub{font-family:var(--mono);font-size:11px;color:var(--tx2);line-height:1.6;display:block;animation:txtIn .36s ease .18s both}

/* ── Viz ── */
.plot-placeholder{display:flex;align-items:center;justify-content:center;flex:1;min-height:400px;color:var(--tx3);font-family:var(--mono);font-size:12px;text-align:center;line-height:2.4}
.viz-panel{animation:vizIn .35s ease both;flex:1;display:flex;flex-direction:column;}
.shiny-plot-output img{animation:plotFade .5s cubic-bezier(.22,1,.36,1);border-radius:var(--r2);display:block}

/* ── Modal ── */
.modal-content{background:var(--srf)!important;border:1px solid var(--bdr2)!important;border-top:2px solid var(--acc2)!important;border-radius:var(--r3)!important;color:var(--tx)!important;box-shadow:var(--sh2)!important}
.modal-header{border-bottom:1px solid var(--bdr)!important;padding:14px 20px!important}
.modal-body{padding:16px 20px!important}
.modal-footer{border-top:1px solid var(--bdr)!important;padding:12px 20px!important;background:var(--srf2)!important;border-radius:0 0 var(--r3) var(--r3)!important;display:flex;gap:8px;justify-content:flex-end}
.modal-title{font-size:14px!important;font-weight:700!important;color:var(--tx)!important}
.close{color:var(--tx3)!important;opacity:.7!important;font-size:18px!important}
.s-section{font-size:10px;font-weight:700;letter-spacing:.14em;text-transform:uppercase;color:var(--tx3);margin:18px 0 10px;padding-top:14px;border-top:1px solid var(--bdr)}
.s-section:first-child{margin-top:2px;border-top:none;padding-top:0}
.s-grid{display:grid;grid-template-columns:1fr 1fr;gap:12px 16px}
.s-field{display:flex;flex-direction:column;gap:5px}
.s-label{font-size:10px;font-weight:700;letter-spacing:.12em;text-transform:uppercase;color:var(--tx2)}
.s-hint{font-size:11px;color:#9ecfb5;line-height:1.5}
.s-info{background:var(--srf2);border:1px solid var(--bdr);border-radius:var(--r2);padding:9px 12px;font-size:11px;color:#9ecfb5;font-family:var(--mono);margin-top:12px;line-height:1.8}
.toggle-row{display:flex;align-items:flex-start;gap:12px;background:var(--srf2);border:1px solid var(--bdr);border-radius:var(--r2);padding:12px 14px;transition:background .15s,border-color .15s}
.toggle-row.on{background:var(--accb);border-color:var(--accm)}
.toggle-desc strong{color:var(--tx);font-size:13px;font-weight:600}
.toggle-desc p{font-size:12px;color:#9ecfb5;margin-top:4px;line-height:1.6}
label{color:var(--tx2)!important}
input[type=number],input[type=text]{background:var(--srf2)!important;color:var(--tx)!important;border:1px solid var(--bdr2)!important;border-radius:var(--r1)!important;font-family:var(--mono)!important;font-size:13px!important;padding:7px 10px!important;transition:border-color .14s!important;width:100%}
input:focus{border-color:var(--acc)!important;outline:none!important;box-shadow:0 0 0 2px rgba(61,153,98,.12)!important}
input[type=checkbox]{accent-color:var(--acc);width:15px;height:15px;margin-top:3px}
.shiny-input-container{margin-bottom:0!important}
select{background:var(--srf2)!important;color:var(--tx)!important;border:1px solid var(--bdr2)!important;border-radius:var(--r1)!important;font-family:var(--mono)!important;font-size:13px!important;padding:7px 10px!important;width:100%}
.fr-val,.risk-pct,.alloc-b{display:inline-block;animation:numPop .3s cubic-bezier(.34,1.56,.64,1) both}
::selection{background:rgba(61,153,98,.35);color:var(--tx)}
button:focus-visible{outline:2px solid var(--acc)!important;outline-offset:2px!important;}

/* ── Results ── */
.results-screen{display:flex;flex-direction:column;gap:20px;padding:4px 0 32px;animation:screenIn .6s cubic-bezier(.22,1,.36,1) both}
.rs-header{text-align:center;padding:36px 20px 24px;background:var(--srf);border:1px solid var(--bdr2);border-top:3px solid var(--acc);border-radius:var(--r3);box-shadow:var(--sh2);}
.rs-title{font-size:34px;font-weight:800;color:var(--tx);letter-spacing:-.01em;margin-bottom:10px;}
.rs-subtitle{font-family:var(--mono);font-size:12px;color:var(--tx3);letter-spacing:.08em;text-transform:uppercase;}
.rs-section-title{font-size:10px;font-weight:700;letter-spacing:.16em;text-transform:uppercase;color:var(--tx3);padding-bottom:8px;border-bottom:1px solid var(--bdr);}
.rs-settings-grid{display:grid;grid-template-columns:repeat(3,1fr);gap:8px;}
.rs-setting{background:var(--srf);border:1px solid var(--bdr);border-radius:var(--r2);padding:10px 14px}
.rs-set-k{font-size:10px;font-weight:700;letter-spacing:.10em;text-transform:uppercase;color:var(--tx3);margin-bottom:3px}
.rs-set-v{font-family:var(--mono);font-size:14px;font-weight:600;color:var(--tx)}
.rs-stat-grid{display:grid;grid-template-columns:repeat(5,1fr);gap:8px}
.rs-stat-card{background:var(--srf);border:1px solid var(--bdr);border-radius:var(--r2);padding:18px 10px;text-align:center;animation:rsStatIn .55s cubic-bezier(.34,1.56,.64,1) both;}
.rs-stat-val{font-family:var(--mono);font-size:28px;font-weight:800;color:var(--tx);line-height:1;margin-bottom:8px;animation:rsNumIn .6s cubic-bezier(.34,1.56,.64,1) .2s both}
.rs-stat-val.green{color:var(--grn)}.rs-stat-val.amber{color:var(--amb)}
.rs-stat-label{font-size:10px;font-weight:600;letter-spacing:.08em;text-transform:uppercase;color:var(--tx3);line-height:1.5;white-space:pre-line}
.rs-trial-list{display:flex;flex-direction:column;gap:6px}
.rs-trial-row{display:grid;grid-template-columns:28px 70px 100px 1fr 1fr 1fr;align-items:center;gap:10px;background:var(--srf);border:1px solid var(--bdr);border-radius:var(--r2);padding:10px 14px;animation:rsRowIn .42s cubic-bezier(.22,1,.36,1) both;}
.rs-trial-row.win{border-left:3px solid var(--grn)}.rs-trial-row.lose{border-left:3px solid var(--red)}.rs-trial-row.det{border-left:3px solid var(--amb)}
.rs-trial-icon{font-size:14px;font-weight:700;text-align:center}
.rs-trial-row.win .rs-trial-icon{color:var(--grn)}.rs-trial-row.lose .rs-trial-icon{color:var(--red)}.rs-trial-row.det .rs-trial-icon{color:var(--amb)}
.rs-trial-num{font-family:var(--mono);font-size:12px;font-weight:600;color:var(--tx)}
.rs-trial-label{font-size:11px;font-weight:600;color:var(--tx2)}
.rs-trial-stat{font-family:var(--mono);font-size:11px;color:var(--tx3)}.rs-trial-stat.err{color:var(--tx2)}
.rs-dl-row{display:grid;grid-template-columns:1fr 1fr;gap:10px}
.rs-dl-card{background:var(--srf);border:1px solid var(--bdr);border-radius:var(--r3);padding:20px 22px;display:flex;flex-direction:column;gap:6px;}
.rs-dl-icon{font-size:22px;margin-bottom:4px}.rs-dl-label{font-size:15px;font-weight:700;color:var(--tx)}.rs-dl-desc{font-size:12px;color:var(--tx2);margin-bottom:8px;line-height:1.5}.rs-dl-btn{margin-top:auto;width:100%!important}
    "))
  ),

  div(class = "wrap",
    div(class = "topbar",
      div(class = "topbar-brand",
        span(class = "topbar-title", "Bayesian Coverage Allocator")
      ),
      div(class = "topbar-right",
        uiOutput("trialPill"),
        uiOutput("graphToggleBtn"),
        uiOutput("settingsBtn"),
        actionButton("rules", "? Rules", class = "topbar-btn")
      )
    ),
    uiOutput("pageContent")
  )
)

# server
server <- function(input, output, session) {

  rv <- reactiveValues(
    screen          = "consent",
    session_id      = paste0("s_", format(Sys.time(), "%Y%m%d_%H%M%S"),
                              "_", sample(10000:99999, 1)),
    game_code       = "",
    is_admin        = FALSE,
    consented       = FALSE,
    viz_mode        = "off",
    hist_type       = "density",  # "density" = old Beta curve bars, "silo" = per-silo bars
    hist_bins       = 20,         # number of bins for density histogram mode
    show_walkthrough = TRUE,      # admin can toggle walkthrough on/off
    show_tutorial    = FALSE,     # admin can toggle pre-practice allocation tutorial
    tutorial_index   = 0,        # which tutorial stage (1,2,3)
    tutorial_done    = FALSE,    # has tutorial been completed
    scripted_mode    = FALSE,    # admin can force fixed hit/miss sequences
    # Settings (admin can change)
    n_silos         = N_SILOS,
    start_units     = START_UNITS,
    total_missiles  = TOTAL_MISSILES,
    risk_k          = RISK_K,
    risk_exp        = RISK_EXP,
    random_missiles = FALSE,
    # Trial state
    scenarios       = make_scenarios(N_TRIALS + 1),
    trial_index     = 0,          # 0 = practice, 1..N_TRIALS = real
    trial_id        = "",
    trial_start_time = Sys.time(),
    a1=1, b1=1, a2=1, b2=1,
    nA=0, kA=0, nB=0, kB=0,
    game_over       = FALSE,
    outcome         = NULL,
    trial_history   = list()
  )

  # derived reactives
  total_s      <- reactive(rv$nA + rv$nB)
  current_risk <- reactive(detection_risk(total_s(), rv$risk_k, rv$risk_exp))
  is_practice  <- reactive(rv$trial_index == 0)
  get_sc       <- function() rv$scenarios[[rv$trial_index + 1]]

  upd_slider <- function() {
    U  <- rv$start_units
    # Keep A's current value if valid, otherwise reset to midpoint
    vA <- if(is.null(input$uA)) floor(U/2) else min(as.integer(input$uA), U)
    vB <- U - vA
    updateSliderInput(session, "uA",        min=0, max=U, value=vA)
    updateSliderInput(session, "uB_slider", min=0, max=U, value=vB)
  }

  # page router
  output$pageContent <- renderUI({
    if      (rv$screen == "consent")  uiOutput("consentScreen")
    else if (rv$screen == "tutorial") uiOutput("tutorialScreen")
    else if (rv$screen == "results")  uiOutput("resultsScreen")
    else {
      tagList(
        div(class="banner-slot", uiOutput("outcomeBanner")),
        div(class="main-grid",
          div(class="left-col",
            div(class="card",
              div(class="card-title","Session Info"),
              uiOutput("settingsDisplay"),
              div(class="divider"),
              uiOutput("budgetDisplay")
            ),
            div(class="card",
              div(class="card-title","Gather Intel"),
              uiOutput("riskGauge"),
              div(class="field-row red",
                span(class="fr-tag red","Field A"),
                span(class="fr-val red",textOutput("pA",inline=TRUE)),
                div(class="fr-meta",
                  div(HTML("<strong>n</strong>&thinsp;"),textOutput("nA",inline=TRUE)),
                  div(HTML("<strong>hits</strong>&thinsp;"),textOutput("kA",inline=TRUE))
                )
              ),
              div(class="field-row blue",
                span(class="fr-tag blue","Field B"),
                span(class="fr-val blue",textOutput("pB",inline=TRUE)),
                div(class="fr-meta",
                  div(HTML("<strong>n</strong>&thinsp;"),textOutput("nB",inline=TRUE)),
                  div(HTML("<strong>hits</strong>&thinsp;"),textOutput("kB",inline=TRUE))
                )
              ),
              div(style="margin-top:14px;",
                actionButton("sampleBoth", "Sample Both Fields", class="btn-sample")
              )
            ),
            div(class="card",
              div(class="card-title","Allocate Coverage"),
              # Field A slider
              div(style="margin-bottom:4px;",
                div(style="display:flex;justify-content:space-between;align-items:baseline;margin-bottom:4px;",
                  span(style="font-size:10px;font-weight:700;letter-spacing:.10em;text-transform:uppercase;color:#ff8096;","Field A"),
                  span(style="font-family:var(--mono);font-size:16px;font-weight:700;color:#ff8096;",textOutput("uA_display",inline=TRUE))
                ),
                div(class="alloc-row",
                  sliderInput("uA",NULL,min=0,max=rv$start_units,
                              value=floor(rv$start_units/2),step=1,width="100%")
                )
              ),
              # Field B slider (mirrors A — always sums to start_units)
              div(style="margin-bottom:12px;",
                div(style="display:flex;justify-content:space-between;align-items:baseline;margin-bottom:4px;",
                  span(style="font-size:10px;font-weight:700;letter-spacing:.10em;text-transform:uppercase;color:#5bc4f0;","Field B"),
                  span(style="font-family:var(--mono);font-size:16px;font-weight:700;color:#5bc4f0;",textOutput("uB_display",inline=TRUE))
                ),
                div(class="alloc-row",
                  sliderInput("uB_slider",NULL,min=0,max=rv$start_units,
                              value=ceiling(rv$start_units/2),step=1,width="100%")
                )
              ),
              actionButton("launch","Launch \u2192",class="btn-primary",width="100%"),
              br(),br(),
              actionButton("nextTrial","Next Trial \u2192",class="btn-next",width="100%")
            )
          ),
          div(class="chart-card",
            div(style="display:flex;justify-content:space-between;align-items:center;margin-bottom:14px;",
              span(class="card-title",style="margin:0;","Posterior Distributions"),
              span(style="font-family:var(--mono);font-size:10px;color:var(--tx3);",
                   "Beta(\u03b1,\u03b2) \u00b7 Bernoulli updates")
            ),
            uiOutput("plotData"),
            uiOutput("plotPanel")
          )
        )
      )
    }
  })

  # consent and code entry screen
  output$consentScreen <- renderUI({
    div(class="consent-screen",
      div(class="consent-card",
        div(class="consent-title","Bayesian Coverage Allocator"),
        div(class="consent-sub","Decision-Making Under Information Constraints"),
        div(class="consent-divider"),

        div(class="consent-section","Game Code"),
        div(class="consent-code-row",
          div(style="font-size:10px;font-weight:700;letter-spacing:.12em;text-transform:uppercase;color:var(--tx2);margin-bottom:4px;",
              "Enter your assigned code"),
          textInput("game_code_input",NULL,value="",placeholder="e.g. A3X7KP"),
          div(class="consent-hint","Your code determines your version of the game. Ask the researcher if you don't have one.")
        ),

        div(class="consent-divider"),

        div(class="consent-section","Data Consent"),
        div(class="consent-body",
          strong("What is collected:"), " how you play — your decisions, sample counts, timing and outcomes.",
          br(),br(),
          strong("How it is used:"), " anonymously, for academic research on decision-making under uncertainty. Your game code is not linked to your identity.",
          br(),br(),
          "You may decline. The game still runs exactly the same — your data will simply not be recorded."
        ),

        div(class="consent-btns",
          actionButton("consentYes","I Consent \u2192 Start",class="btn-primary"),
          actionButton("consentNo", "Decline \u2014 Play Anyway",class="btn-decline")
        )
      )
    )
  })

  observeEvent(input$consentYes, { handle_code(consented=TRUE)  })
  observeEvent(input$consentNo,  { handle_code(consented=FALSE) })

  handle_code <- function(consented) {
    raw  <- input$game_code_input %||% ""
    code <- toupper(trimws(raw))

    if (nchar(code) == 0) {
      showNotification("Please enter your game code first.", type="warning", duration=3)
      return()
    }

    # admin bypass
    if (tolower(code) == ADMIN_CODE) {
      rv$is_admin  <- TRUE
      rv$consented <- consented
      rv$game_code <- ADMIN_CODE
      rv$viz_mode  <- "curve"
      show_admin_modal()
      return()
    }

    # look up viz condition from code
    viz <- lookup_code(code)
    if (is.null(viz)) {
      showNotification(paste0("Code '",code,"' not recognised. Please check with the researcher."),
                       type="error", duration=4)
      return()
    }

    rv$game_code       <- code
    rv$consented       <- consented
    rv$is_admin        <- FALSE
    rv$viz_mode        <- viz
    rv$hist_type       <- "density"
    rv$show_walkthrough <- TRUE
    rv$show_tutorial    <- TRUE
    rv$scripted_mode    <- FALSE
    rv$scenarios       <- make_scenarios(N_TRIALS + 1)

    rv$tutorial_index <- 1
    rv$tutorial_done  <- FALSE
    rv$screen         <- "tutorial"
  }

  # admin modal
  show_admin_modal <- function() {
    showModal(modalDialog(
      title="Admin Settings", easyClose=FALSE, footer=NULL,
      div(
        div(class="s-section","Game Parameters"),
        div(class="s-grid",
          div(class="s-field",div(class="s-label","Coverage Units"),
              numericInput("adm_units",NULL,value=rv$start_units,min=1,max=NA,step=1,width="100%"),
              div(class="s-hint","Must be >= 1")),
          div(class="s-field",div(class="s-label","Total Missiles"),
              numericInput("adm_missiles",NULL,value=rv$total_missiles,min=1,max=NA,step=1,width="100%"),
              div(class="s-hint","Must be < silos per field")),
          div(class="s-field",div(class="s-label","Silos per Field"),
              numericInput("adm_silos",NULL,value=rv$n_silos,min=1,max=NA,step=1,width="100%"),
              div(class="s-hint","Must be > total missiles")),
          div(class="s-field",div(class="s-label","Visualisation"),
              selectInput("adm_viz",NULL,
                choices=c("Off"="off","Curve"="curve","Interval Bar"="dotplot",
                          "Gauge"="gauge","Histogram"="histogram"),
                selected=rv$viz_mode,width="100%"))
        ),
        div(class="s-section","Histogram Settings (when Histogram viz is selected)"),
        div(class="s-grid",
          div(class="s-field",div(class="s-label","Histogram Type"),
              selectInput("adm_hist_type",NULL,
                choices=c("Per-silo bars (new)"="silo","Density bars (original)"="density"),
                selected=rv$hist_type,width="100%"),
              div(class="s-hint","Silo: one bar per silo. Density: Beta distribution shape.")),
          div(class="s-field",div(class="s-label","Density Bars (density mode only)"),
              numericInput("adm_hist_bins",NULL,value=rv$hist_bins,min=5,max=NA,step=5,width="100%"),
              div(class="s-hint","Number of bars in density histogram."))
        ),
        div(class="s-section","Risk Curve"),
        div(class="s-grid",
          div(class="s-field",div(class="s-label","Steepness k"),
              numericInput("adm_k",NULL,value=rv$risk_k,min=0.0001,max=NA,step=0.005,width="100%"),
              div(class="s-hint","Default: 0.001")),
          div(class="s-field",div(class="s-label","Exponent"),
              numericInput("adm_exp",NULL,value=rv$risk_exp,min=0.1,max=NA,step=0.1,width="100%"),
              div(class="s-hint","Default: 2.5"))
        ),
        div(class="s-section","Study Experience"),
        div(class=if(rv$show_walkthrough)"toggle-row on"else"toggle-row",
          checkboxInput("adm_walkthrough",NULL,value=rv$show_walkthrough),
          div(class="toggle-desc",
            strong("Show Walkthrough"),
            p("On: participants see a step-by-step tour when the game starts. Off: skip straight to play.")
          )
        ),
        div(style="margin-top:8px;",
          div(class=if(rv$show_tutorial)"toggle-row on"else"toggle-row",
            checkboxInput("adm_tutorial",NULL,value=rv$show_tutorial),
            div(class="toggle-desc",
              strong("Allocation Tutorial"),
              p("On: 3 guided rounds before practice where probabilities are revealed and the player just allocates — no sampling required.")
            )
          )
        ),
        div(style="margin-top:8px;",
          div(class=if(rv$scripted_mode)"toggle-row on"else"toggle-row",
            checkboxInput("adm_scripted",NULL,value=rv$scripted_mode),
            div(class="toggle-desc",
              strong("Scripted Graph Mode"),
              p("On: all participants see identical hit/miss sequences — graphs change the same way every session. Off: sequences are randomly sampled from true probabilities.")
            )
          )
        ),
        div(class="s-section","Missile Distribution"),
        div(class=if(rv$random_missiles)"toggle-row on"else"toggle-row",
          checkboxInput("adm_random",NULL,value=rv$random_missiles),
          div(class="toggle-desc",
            strong("Randomise Missile Split"),
            p("Off: fixed seeded split (reproducible across all sessions).","On: random each run.")
          )
        ),
        div(class="s-info",paste0(
          "Always risk mode  \u00b7  ",N_TRIALS," real trials + 1 practice  \u00b7  Max ",MAX_SAMPLES," combined samples")),
        br(),
        div(style="display:flex;gap:8px;justify-content:flex-end;",
          actionButton("admCancel","Cancel",class="btn-cancel"),
          actionButton("admStart","Start Study \u2192",class="btn-primary")
        )
      )
    ))
  }

  observeEvent(input$admCancel, { removeModal() })

  observeEvent(input$admStart, {
    units    <- max(1, as.integer(input$adm_units    %||% START_UNITS))
    missiles <- max(1, as.integer(input$adm_missiles %||% TOTAL_MISSILES))
    silos    <- max(1, as.integer(input$adm_silos    %||% N_SILOS))
    k        <- max(0.0001, as.numeric(input$adm_k   %||% RISK_K))
    exp_val  <- max(0.1,    as.numeric(input$adm_exp %||% RISK_EXP))

    # Validate: missiles must be strictly less than silos
    if (missiles >= silos) {
      showNotification(
        paste0("Total missiles (", missiles, ") must be less than silos per field (", silos, "). Please fix before starting."),
        type = "error", duration = 5)
      return()
    }

    rv$start_units     <- units
    rv$total_missiles  <- missiles
    rv$n_silos         <- silos
    rv$viz_mode        <- input$adm_viz %||% "off"
    rv$hist_type       <- input$adm_hist_type %||% "silo"
    rv$hist_bins       <- max(5, as.integer(input$adm_hist_bins %||% 20))
    rv$show_walkthrough <- isTRUE(input$adm_walkthrough)
    rv$show_tutorial    <- isTRUE(input$adm_tutorial)
    rv$scripted_mode    <- isTRUE(input$adm_scripted)
    rv$risk_k          <- k
    rv$risk_exp        <- exp_val
    rv$random_missiles <- isTRUE(input$adm_random)
    rv$scenarios       <- make_scenarios(N_TRIALS + 1, rv$total_missiles)
    removeModal()
    if (rv$show_tutorial) {
      rv$tutorial_index <- 1
      rv$tutorial_done  <- FALSE
      rv$screen <- "tutorial"
    } else {
      rv$screen <- "game"
    # Inject step-by-step walkthrough after a short delay so game UI renders first
    start_trial()
    }  # end else (no tutorial)
  })

  # viz toggle button for admin
  output$graphToggleBtn <- renderUI({
    if (!rv$is_admin) return(NULL)
    lbl <- switch(rv$viz_mode,
      "off"="Viz: Off","curve"="Viz: Curve","dotplot"="Viz: Interval",
      "gauge"="Viz: Gauge","histogram"="Viz: Histogram","Viz: Off")
    cls <- if(rv$viz_mode=="off")"topbar-btn"else"topbar-btn on"
    actionButton("toggleGraph", lbl, class=cls)
  })
  observeEvent(input$toggleGraph, {
    if (!rv$is_admin) return()
    rv$viz_mode <- switch(rv$viz_mode,
      "off"="curve","curve"="dotplot","dotplot"="gauge",
      "gauge"="histogram","histogram"="off","off")
  })

  # settings button for admin
  output$settingsBtn <- renderUI({
    if (!rv$is_admin) return(NULL)
    actionButton("openSettings","Settings",class="topbar-btn")
  })
  observeEvent(input$openSettings, { if (rv$is_admin) show_admin_modal() })

  # trial progress pill in topbar
  output$trialPill <- renderUI({
    if (rv$screen != "game") return(span(class="pill","Ready"))
    if (is_practice()) span(class="pill practice-pill","Practice Trial")
    else span(class="pill active",paste0("Trial ",rv$trial_index," / ",N_TRIALS))
  })

  # trial setup and teardown
  fire_walkthrough <- function() {
    shinyjs::delay(700, {
      shinyjs::runjs("
(function(){
  if(window.WT) return; // already running
  var steps=[
    {dim:[],
     title:'Welcome',
     text:'Two fields hide missiles. Identify the more dangerous one, then deploy your units.',
    },
    {dim:['.left-col .card:nth-child(1)','.left-col .card:nth-child(3)','.chart-card','.topbar','.banner-slot'],
     title:'Field Intel',
     text:'Each row shows your probability estimate for that field. Both start at 0.500 — equal uncertainty.',
    },
    {dim:['.left-col .card:nth-child(1)','.left-col .card:nth-child(3)','.chart-card','.topbar','.banner-slot','.field-row.red','.field-row.blue'],
     title:'Sample Both Fields',
     text:'Press this to take one intel reading from each field. Up to 8 presses per round.',
    },
    {dim:['.left-col .card:nth-child(1)','.left-col .card:nth-child(3)','.chart-card','.topbar','.banner-slot','.field-row.red','.field-row.blue'],
     title:'Detonation Risk',
     text:'Every sample raises this. First press is always safe. Gets serious around press 6.',
    },
    {dim:['.left-col','.topbar','.banner-slot'],
     title:'Visualization',
     text:'This updates after every sample. Use it to judge which field needs more coverage.',
    },
    {dim:['.left-col .card:nth-child(1)','.left-col .card:nth-child(2)','.chart-card','.topbar','.banner-slot'],
     title:'Allocate',
     text:'Split your units between the two fields using these sliders. They always sum to your total.',
    },
    {dim:['.left-col .card:nth-child(1)','.left-col .card:nth-child(2)','.chart-card','.topbar','.banner-slot'],
     title:'Launch',
     text:'When ready, hit Launch. Cover all missiles in both fields to win the round.',
    },
    {dim:[],
     title:'Practice Round',
     text:'This round is practice — nothing is recorded. Get comfortable, then the real trials begin.',
    }
  ];
  var st=document.createElement('style');
  st.textContent=[
    '.wt-dim{transition:filter .35s ease!important;filter:brightness(0.22) saturate(0.4)!important;}',
    '@keyframes wtCardIn{from{opacity:0;transform:translateY(10px)}to{opacity:1;transform:none}}',
    '@keyframes wtCardOut{from{opacity:1;transform:none}to{opacity:0;transform:translateY(-8px)}}'
  ].join('');
  document.head.appendChild(st);
  var card=document.createElement('div');
  card.style.cssText=[
    'position:fixed;bottom:32px;left:50%;transform:translateX(-50%)',
    'z-index:9999;width:340px',
    'background:#0c1c10;border:1px solid #2d5c3a;border-top:3px solid #5ecf8a',
    'border-radius:16px;padding:24px 28px 20px',
    'box-shadow:0 20px 60px rgba(0,0,0,.90)',
    'pointer-events:auto',
    'animation:wtCardIn .38s cubic-bezier(.22,1,.36,1) both'
  ].join(';');
  document.body.appendChild(card);
  function render(idx){
    var s=steps[idx];
    var isFirst=idx===0,isLast=idx===steps.length-1;
    document.querySelectorAll('.wt-dim').forEach(function(e){e.classList.remove('wt-dim');});
    s.dim.forEach(function(sel){document.querySelectorAll(sel).forEach(function(e){e.classList.add('wt-dim');});});
    var dots='';
    for(var i=0;i<steps.length;i++)
      dots+='<span style=\"display:inline-block;width:6px;height:6px;border-radius:50%;margin:0 3px;background:'+(i===idx?'#5ecf8a':'#1e3828')+';transition:background .25s\"></span>';
    var back=isFirst?'':'<button onclick=\"WT.go('+(idx-1)+')\" style=\"flex:1;background:transparent;border:1px solid #2d5c3a;color:#7ec49a;border-radius:9px;padding:9px 0;font-size:12px;font-weight:600;cursor:pointer;font-family:Inter,sans-serif\">← Back</button>';
    var next='<button onclick=\"'+(isLast?'WT.done()':'WT.go('+(idx+1)+')')+'\" style=\"flex:2;background:#3d9962;border:none;color:#030a05;border-radius:9px;padding:9px 0;font-size:13px;font-weight:700;cursor:pointer;font-family:Inter,sans-serif\">'+(isLast?'Start Playing':'Next')+'</button>';
    card.innerHTML=[
      '<div style=\"display:flex;align-items:center;gap:12px;margin-bottom:10px\">',
        '<span style=\"font-size:14px;font-weight:800;color:#e0f2e8;font-family:Inter,sans-serif\">'+s.title+'</span>',
        '<span style=\"margin-left:auto;font-size:11px;color:#3d6b50;font-family:Inter,sans-serif\">'+(idx+1)+'/'+steps.length+'</span>',
      '</div>',
      '<div style=\"font-size:12.5px;color:#9ecfb5;line-height:1.65;margin-bottom:16px;font-family:Inter,sans-serif\">'+s.text+'</div>',
      '<div style=\"text-align:center;margin-bottom:14px\">'+dots+'</div>',
      '<div style=\"display:flex;gap:8px\">'+back+next+'</div>'
    ].join('');
  }
  // block all app interaction while walkthrough is active
  var blocker=document.createElement('div');
  blocker.id='wt-blocker';
  blocker.style.cssText='position:fixed;inset:0;z-index:9998;pointer-events:auto;';
  document.body.appendChild(blocker);

  window.WT={
    go:function(idx){render(idx);},
    done:function(){
      document.querySelectorAll('.wt-dim').forEach(function(e){e.classList.remove('wt-dim');});
      card.style.animation='wtCardOut .3s cubic-bezier(.22,1,.36,1) forwards';
      var bl=document.getElementById('wt-blocker');
      if(bl) bl.remove();
      setTimeout(function(){card.remove();st.remove();delete window.WT;},320);
    }
  };
  render(0);
})();
      ")
    })
  }

  start_trial <- function() {
    rv$trial_id         <- paste0(rv$session_id,"_t",rv$trial_index)
    rv$trial_start_time <- Sys.time()
    rv$a1<-1;rv$b1<-1;rv$a2<-1;rv$b2<-1
    rv$nA<-0;rv$kA<-0;rv$nB<-0;rv$kB<-0
    rv$game_over<-FALSE; rv$outcome<-NULL
    updateSliderInput(session,"uA",       min=0,max=rv$start_units,value=floor(rv$start_units/2))
    updateSliderInput(session,"uB_slider",min=0,max=rv$start_units,value=ceiling(rv$start_units/2))
    shinyjs::disable("nextTrial")
    shinyjs::enable("sampleBoth"); shinyjs::enable("launch")
    if (rv$consented && !is_practice()) log_ev("trial_start")
    # show walkthrough at the start of the practice trial
    if (rv$trial_index == 0 && rv$show_walkthrough) fire_walkthrough()
  }

  trial_complete <- function(player_uA = NULL) {
    shinyjs::disable("sampleBoth"); shinyjs::disable("launch")
    o <- rv$outcome
    if (!is_practice()) {
      p_uA <- if (!is.null(player_uA)) player_uA else o$uA
      U    <- rv$start_units

      # true optimal = exact missile counts, or proportional if budget is tight
      mA <- o$mA; mB <- o$mB
      if (!is.null(mA) && !is.null(mB) && !is.na(mA) && !is.na(mB)) {
        if (mA + mB <= U) {
          true_opt_A <- mA
          true_opt_B <- mB
        } else {
          # round and give remainder to heavier field
          true_opt_A <- round(U * mA / (mA + mB))
          true_opt_B <- U - true_opt_A
        }
        true_split_err <- abs(p_uA - true_opt_A)
      } else {
        true_opt_A <- NA; true_opt_B <- NA; true_split_err <- NA
      }

      rv$trial_history <- c(rv$trial_history, list(list(
        trial          = rv$trial_index,
        win            = isTRUE(o$win),
        detected       = isTRUE(o$detected),
        samples        = total_s(),
        missed         = o$missed_total,
        player_uA      = p_uA,
        player_uB      = U - p_uA,
        true_opt_A     = true_opt_A,
        true_opt_B     = true_opt_B,
        split_err      = if(!isTRUE(o$detected)) true_split_err else NA,
        elapsed        = round(as.numeric(difftime(Sys.time(),rv$trial_start_time,units="secs")),1)
      )))
    }
    if (rv$trial_index < N_TRIALS) shinyjs::enable("nextTrial")
    else { shinyjs::disable("nextTrial"); rv$screen <- "results" }
  }

  # sample button handler
  # ONE detonation roll per button press, using the count BEFORE this press.
  # This means press 1 (s=0 before) always has 0% risk.
  # Both fields are sampled regardless — the risk is for the whole press.
  observeEvent(input$sampleBoth, {
    req(!rv$game_over, rv$nA < MAX_SAMPLES)
    sc <- get_sc()

    # get risk before sampling so press 1 is always 0%
    pre_risk <- detection_risk(total_s(), rv$risk_k, rv$risk_exp)

    # pick hit/miss from scripted or random stream
    stream_idx <- min(rv$trial_index + 1, length(SCRIPTED_STREAMS))
    dA <- if (rv$scripted_mode)
      SCRIPTED_STREAMS[[stream_idx]]$A[rv$nA + 1]
    else
      sc$sample_A_stream[rv$nA + 1]

    dB <- if (rv$scripted_mode)
      SCRIPTED_STREAMS[[stream_idx]]$B[rv$nB + 1]
    else
      sc$sample_B_stream[rv$nB + 1]

    rv$a1 <- rv$a1 + dA; rv$b1 <- rv$b1 + (1-dA); rv$nA <- rv$nA+1; rv$kA <- rv$kA+dA
    rv$a2 <- rv$a2 + dB; rv$b2 <- rv$b2 + (1-dB); rv$nB <- rv$nB+1; rv$kB <- rv$kB+dB

    upd_slider()
    if (rv$consented && !is_practice()) log_ev("sample_both", dA, dB)

    # Single roll after both samples
    if (runif(1) < pre_risk) end_detonation(pre_risk)
  })

  end_detonation <- function(pre_risk=NULL) {
    sc  <- get_sc()
    pAp <- rv$a1/(rv$a1+rv$b1); pBp <- rv$a2/(rv$a2+rv$b2)
    U   <- rv$start_units; uA <- min(as.integer(isolate(input$uA)), U)
    sr  <- if(!is.null(pre_risk)) round(pre_risk*100,1) else round(current_risk()*100,1)
    rv$game_over <- TRUE
    rv$outcome   <- list(
      detected=TRUE, win=FALSE, risk_pct=sr, U=U, uA=uA, uB=U-uA,
      mA=sc$mA_true, mB=sc$mB_true, hitA=0, hitB=0,
      missed_A=sc$mA_true, missed_B=sc$mB_true, missed_total=sc$mA_true+sc$mB_true,
      pA_true=sc$pA_true, pB_true=sc$pB_true,
      pA_post=round(pAp,3), pB_post=round(pBp,3), opt_uA=NA, opt_uB=NA)
    if (rv$consented && !is_practice()) {
      r <- build_launch_row()
      r$detected     <- TRUE
      r$success      <- FALSE
      r$missed_a     <- sc$mA_true
      r$missed_b     <- sc$mB_true
      r$missed_total <- sc$mA_true + sc$mB_true
      log_row(r,"launch_trials",LAUNCHES_FILE)
    }
    trial_complete()
  }

  observeEvent(input$launch, {
    req(!rv$game_over)
    sc  <- get_sc(); U <- rv$start_units
    uA  <- min(as.integer(input$uA),U); uB <- U-uA
    pAp <- rv$a1/(rv$a1+rv$b1); pBp <- rv$a2/(rv$a2+rv$b2)
    hitA<-min(sc$mA_true,uA);hitB<-min(sc$mB_true,uB)
    mA  <-max(sc$mA_true-uA,0);mB<-max(sc$mB_true-uB,0)
    opt <- bayes_optimal_split(rv$a1,rv$b1,rv$a2,rv$b2,rv$n_silos,U)
    rv$game_over<-TRUE
    rv$outcome<-list(
      detected=FALSE, win=(mA+mB)==0, risk_pct=round(current_risk()*100,1),
      U=U, uA=uA, uB=uB, mA=sc$mA_true, mB=sc$mB_true,
      hitA=hitA, hitB=hitB, missed_A=mA, missed_B=mB, missed_total=mA+mB,
      pA_true=sc$pA_true, pB_true=sc$pB_true,
      pA_post=round(pAp,3), pB_post=round(pBp,3),
      opt_uA=opt$opt_uA, opt_uB=opt$opt_uB)
    if (rv$consented && !is_practice()) {
      r <- build_launch_row(opt=opt)
      # All column names must be lowercase to match Supabase table
      r$detected     <- FALSE
      r$success      <- rv$outcome$win
      r$alloc_a      <- uA
      r$alloc_b      <- uB
      r$missed_a     <- mA
      r$missed_b     <- mB
      r$missed_total <- mA + mB
      log_row(r,"launch_trials",LAUNCHES_FILE)
    }
    trial_complete(player_uA = uA)
  })

  observeEvent(input$nextTrial, {
    req(rv$game_over, rv$trial_index <= N_TRIALS)
    rv$trial_index <- rv$trial_index + 1
    start_trial()
  })

  observe({
    if (rv$game_over || rv$nA >= MAX_SAMPLES) shinyjs::disable("sampleBoth")
    else shinyjs::enable("sampleBoth")
    if (rv$game_over) shinyjs::disable("launch") else shinyjs::enable("launch")
  })

  # build rows for logging
  build_base_row <- function() {
    sc  <- get_sc()
    pAp <- rv$a1/(rv$a1+rv$b1); pBp <- rv$a2/(rv$a2+rv$b2)
    list(
      timestamp      = as.character(Sys.time()),
      session_id     = rv$session_id,
      game_code      = rv$game_code,
      viz_mode       = rv$viz_mode,
      trial_id       = rv$trial_id,
      trial_index    = rv$trial_index,
      is_practice    = is_practice(),
      n_silos        = rv$n_silos,
      start_units    = rv$start_units,
      total_missiles = rv$total_missiles,
      risk_k         = rv$risk_k,
      risk_exp       = rv$risk_exp,
      pa_true        = sc$pA_true,
      pb_true        = sc$pB_true,
      ma_true        = sc$mA_true,
      mb_true        = sc$mB_true,
      a1=rv$a1, b1=rv$b1, a2=rv$a2, b2=rv$b2,
      post_mean_a    = round(pAp,4),
      post_mean_b    = round(pBp,4),
      post_var_total = round(posterior_var(rv$a1,rv$b1)+posterior_var(rv$a2,rv$b2),6),
      samp_total     = total_s(),
      detection_risk = round(current_risk()*100,2),
      elapsed        = round(as.numeric(difftime(Sys.time(),rv$trial_start_time,units="secs")),2)
    )
  }

  build_launch_row <- function(opt=NULL) {
    r <- build_base_row()
    r$event_type    <- "launch"
    r$alloc_a       <- min(as.integer(isolate(input$uA)), rv$start_units)
    r$alloc_b       <- rv$start_units - r$alloc_a
    r$opt_a         <- if(!is.null(opt) && !is.na(opt$opt_uA)) opt$opt_uA else NA
    r$opt_b         <- if(!is.null(opt) && !is.na(opt$opt_uB)) opt$opt_uB else NA
    r$alloc_err_abs <- if(!is.null(opt) && !is.na(opt$opt_uA)) abs(r$alloc_a - opt$opt_uA) else NA
    r$h1_samples    <- total_s()
    r$h4_post_var   <- round(posterior_var(rv$a1,rv$b1)+posterior_var(rv$a2,rv$b2),6)
    r
  }

  log_ev <- function(et, sr_a=NA, sr_b=NA) {
    r <- build_base_row()
    r$event_type      <- et
    r$sample_result_a <- sr_a
    r$sample_result_b <- sr_b
    log_row(r,"all_events",EVENTS_FILE)
  }

  # reactive display outputs
  output$settingsDisplay <- renderUI({
    vl <- switch(rv$viz_mode,"off"="None","curve"="Curve","dotplot"="Interval",
                 "gauge"="Gauge","histogram"="Histogram","None")
    div(class="settings-display",HTML(paste0(
      "<span class='badge-risk'>Risk Mode</span>",
      " \u00b7 Units: <strong>",rv$start_units,"</strong>",
      " \u00b7 Missiles: <strong>",rv$total_missiles,"</strong>",
      "<br>Viz: <strong>",vl,"</strong>",
      " \u00b7 k=<strong>",rv$risk_k,"</strong>",
      " \u00b7 exp=<strong>",rv$risk_exp,"</strong>"
    )))
  })

  output$budgetDisplay <- renderUI({
    div(class="status-strip",
      div(class="kv",span(class="k","samples"),span(class="v",total_s())),
      div(class="sep","\u00b7"),
      div(class="kv",span(class="k","units"),span(class="v",rv$start_units)),
      div(class="sep","\u00b7"),
      div(class="kv",span(class="k","risk"),
          span(class="v",style="color:var(--amb);",paste0(round(current_risk()*100,1),"%")))
    )
  })

  output$riskGauge <- renderUI({
    pct  <- current_risk()*100
    pfmt <- sprintf("%.1f%%",pct)
    col  <- if(pct<10)"var(--tx2)"else if(pct<35)"var(--amb)"else"var(--red)"
    cls  <- if(pct>=45)"risk-gauge hi"else"risk-gauge"
    prev <- paste(sapply(1:6,function(s)
      sprintf("s%d->%.0f%%",s,detection_risk(s,rv$risk_k,rv$risk_exp)*100)),collapse="  ")
    div(class=cls,
      div(class="risk-top",
        span(class="risk-title","Detonation Risk"),
        span(class="risk-pct",style=paste0("color:",col,";"),pfmt)
      ),
      div(class="risk-track",
        div(class="risk-fill",style=paste0("width:",min(pct,100),"%;background:",col,";"))),
      div(class="risk-hint","Roll fires after each combined sample.", br(), prev)
    )
  })

  output$plotData <- renderUI({
    tags$div(id="plotData",style="display:none;",
      `data-a1`=rv$a1,`data-b1`=rv$b1,`data-a2`=rv$a2,`data-b2`=rv$b2)
  })

  output$pA <- renderText(sprintf("%.3f",rv$a1/(rv$a1+rv$b1)))
  output$pB <- renderText(sprintf("%.3f",rv$a2/(rv$a2+rv$b2)))
  output$nA <- renderText(paste0(rv$nA,"/",MAX_SAMPLES))
  output$kA <- renderText(rv$kA)
  output$nB <- renderText(paste0(rv$nB,"/",MAX_SAMPLES))
  output$kB <- renderText(rv$kB)

  # Dual slider sync — moving A adjusts B and vice versa
  # Each slider is the authoritative source when it last changed
  observeEvent(input$uA, {
    req(!is.null(input$uA))
    new_uB <- rv$start_units - as.integer(input$uA)
    if (!is.null(input$uB_slider) && as.integer(input$uB_slider) != new_uB)
      updateSliderInput(session,"uB_slider",value=new_uB)
  }, ignoreInit=TRUE)

  observeEvent(input$uB_slider, {
    req(!is.null(input$uB_slider))
    new_uA <- rv$start_units - as.integer(input$uB_slider)
    if (!is.null(input$uA) && as.integer(input$uA) != new_uA)
      updateSliderInput(session,"uA",value=new_uA)
  }, ignoreInit=TRUE)

  output$uA_display <- renderText(as.character(min(input$uA %||% floor(rv$start_units/2), rv$start_units)))
  output$uB_display <- renderText(as.character(rv$start_units - min(input$uA %||% floor(rv$start_units/2), rv$start_units)))
  output$uB <- renderText(paste0("B: ", rv$start_units - min(input$uA %||% 0, rv$start_units)))

  # outcome popup
  output$outcomeBanner <- renderUI({ NULL })   # strip kept empty; modal handles display

  observe({
    if (is.null(rv$outcome)) return()
    o    <- rv$outcome
    prac <- is_practice()
    sub_txt <- if (prac) "<p style='color:#4a8c68;font-size:12px;margin-top:8px;font-family:monospace;'>Practice round — not recorded.</p>" else ""

    if (isTRUE(o$detected)) {
      icon  <- '<span style="font-size:11px;font-weight:800;letter-spacing:.14em;color:#b88a38;font-family:Inter,sans-serif;">DETONATED</span>'
      title_txt <- if(prac) "Practice \u2014 Missile Detonated" else "Missile Detonated!"
      body  <- paste0(
        "<div style='text-align:center;padding:8px 0;'>",
        "<div style='font-size:48px;margin-bottom:12px;'>", icon, "</div>",
        "<div style='font-size:15px;font-weight:600;color:#d8efe2;margin-bottom:6px;'>",
        sprintf("Risk was %.1f%% when the roll fired.", o$risk_pct), "</div>",
        "<div style='font-size:12px;color:#7ec49a;font-family:monospace;'>",
        sprintf("Samples: %d", total_s()), "</div>",
        sub_txt, "</div>")
      border_col <- "#b88a38"

    } else if (o$win) {
      icon  <- '<span style="font-size:11px;font-weight:800;letter-spacing:.14em;color:#3db87a;font-family:Inter,sans-serif;">MISSION SUCCESS</span>'
      title_txt <- if(prac) "Practice \u2014 Mission Success!" else "Mission Success!"
      body  <- paste0(
        "<div style='text-align:center;padding:8px 0;'>",
        "<div style='font-size:48px;margin-bottom:12px;'>", icon, "</div>",
        "<div style='font-size:15px;font-weight:600;color:#d8efe2;margin-bottom:6px;'>All missiles covered!</div>",
        "<div style='font-size:12px;color:#7ec49a;font-family:monospace;'>",
        sprintf("Budget %d  \u00b7  A: %d/%d  \u00b7  B: %d/%d", o$U, o$hitA, o$mA, o$hitB, o$mB), "</div>",
        if (!prac && !is.na(o$opt_uA)) paste0(
          "<div style='font-size:11px;color:#4a8c68;font-family:monospace;margin-top:4px;'>",
          sprintf("Your split A%d/B%d  \u00b7  Optimal A%d/B%d", o$uA, o$uB, o$opt_uA, o$opt_uB),
          "</div>") else "",
        sub_txt, "</div>")
      border_col <- "#3db87a"

    } else {
      icon  <- '<span style="font-size:11px;font-weight:800;letter-spacing:.14em;color:#c96070;font-family:Inter,sans-serif;">MISSION FAILED</span>'
      title_txt <- if(prac) paste0("Practice \u2014 ", o$missed_total, " Missile(s) Escaped")
                   else      paste0(o$missed_total, " Missile(s) Escaped!")
      body  <- paste0(
        "<div style='text-align:center;padding:8px 0;'>",
        "<div style='font-size:48px;margin-bottom:12px;'>", icon, "</div>",
        "<div style='font-size:15px;font-weight:600;color:#d8efe2;margin-bottom:6px;'>",
        sprintf("%d missile(s) got through.", o$missed_total), "</div>",
        "<div style='font-size:12px;color:#7ec49a;font-family:monospace;'>",
        sprintf("Budget %d  \u00b7  A: %d/%d (%d missed)  \u00b7  B: %d/%d (%d missed)",
                o$U, o$hitA, o$mA, o$missed_A, o$hitB, o$mB, o$missed_B), "</div>",
        if (!prac && !is.na(o$opt_uA)) paste0(
          "<div style='font-size:11px;color:#4a8c68;font-family:monospace;margin-top:4px;'>",
          sprintf("Your split A%d/B%d  \u00b7  Optimal A%d/B%d", o$uA, o$uB, o$opt_uA, o$opt_uB),
          "</div>") else "",
        sub_txt, "</div>")
      border_col <- "#c96070"
    }

    showModal(modalDialog(
      title = title_txt,
      easyClose = TRUE,
      footer = modalButton("Continue"),
      tags$div(
        style = paste0("border-top:3px solid ", border_col, ";margin:-1px -1px 0;border-radius:var(--r3) var(--r3) 0 0;"),
        HTML(body)
      )
    ))
  }) |> bindEvent(rv$outcome, ignoreNULL = TRUE, ignoreInit = TRUE)

  # visualization panel
  output$plotPanel <- renderUI({
    mode <- rv$viz_mode

    if (mode == "off") {
      div(class="plot-placeholder",div(
        div(style="color:var(--tx2);font-size:13px;","Visualization hidden"),
        div(style="font-size:11px;color:var(--tx3);margin-top:5px;",
            "Use the numerical estimates to inform your allocation.")
      ))

    } else if (mode == "curve") {
      plotOutput("plot",height="440px")

    } else if (mode == "dotplot") {
      a1<-rv$a1;b1<-rv$b1;a2<-rv$a2;b2<-rv$b2;nA<-rv$nA;nB<-rv$nB
      make_iv_svg <- function(a,b,col,uid){
        mn<-a/(a+b);lo<-qbeta(.05,a,b);hi<-qbeta(.95,a,b)
        sw<-340;sh<-80;pl<-24;pr<-24;pw<-sw-pl-pr
        by<-22;bh<-22;ly<-56;ay<-64;ty<-77
        v2x<-function(v)pl+v*pw
        lx<-v2x(lo);hx<-v2x(hi);mx<-v2x(mn)
        llx<-max(pl+4,round(lx,1));hlx<-min(sw-pr-4,round(hx,1))
        if((hlx-llx)<34){mid<-(llx+hlx)/2;llx<-mid-17;hlx<-mid+17}
        grad_id<-paste0("ivg",uid)
        defs<-tags$defs(
          tags$linearGradient(id=grad_id,x1="0",y1="0",x2="1",y2="0",
            tags$stop(offset="0%",  style=paste0("stop-color:",col,";stop-opacity:0.55")),
            tags$stop(offset="50%", style=paste0("stop-color:",col,";stop-opacity:0.95")),
            tags$stop(offset="100%",style=paste0("stop-color:",col,";stop-opacity:0.55"))
          )
        )
        tags$svg(xmlns="http://www.w3.org/2000/svg",viewBox=paste0("0 0 ",sw," ",sh),style="width:100%;display:block;",
          defs,
          # Track
          tags$rect(x=pl,y=by,width=pw,height=bh,rx="5",fill="#050d08",stroke="#1a3024",`stroke-width`="1"),
          # Interval bar with gradient
          tags$rect(x=round(lx,1),y=by,width=max(5,round(hx-lx,1)),height=bh,rx="4",
                    fill=paste0("url(#",grad_id,")"),opacity="0.88"),
          # White centre marker
          tags$rect(x=round(mx-3,1),y=by-4,width="6",height=bh+8,rx="3",
                    fill="#ffffff",opacity="0.95"),
          # Glow on centre marker
          tags$rect(x=round(mx-6,1),y=by-1,width="12",height=bh+2,rx="6",
                    fill=col,opacity="0.15"),
          # Mean label — 11px, above bar, clear of bar top (by=22, label at y=14)
          tags$text(x=round(mx,1),y=14,`text-anchor`="middle",
            style=paste0("font-family:monospace;font-size:11px;fill:",col,";font-weight:700;"),
            sprintf("%.3f",mn)),
          # Lo/hi boundary labels — 10px at ly=56 (well between bar bottom at 44 and axis at 64)
          tags$text(x=llx,y=ly,`text-anchor`="middle",
            style="font-family:monospace;font-size:10px;fill:#4a8c68;",sprintf("%.2f",lo)),
          tags$text(x=hlx,y=ly,`text-anchor`="middle",
            style="font-family:monospace;font-size:10px;fill:#4a8c68;",sprintf("%.2f",hi)),
          # Axis
          tags$line(x1=pl,y1=ay,x2=sw-pr,y2=ay,stroke="#1a3024",`stroke-width`="1.2"),
          # Tick marks + labels at 77px (clear of axis at 64)
          tagList(lapply(c(0,.25,.5,.75,1),function(v){tx<-v2x(v);tagList(
            tags$line(x1=round(tx,1),y1=ay,x2=round(tx,1),y2=ay+4,stroke="#22402e",`stroke-width`="1"),
            tags$text(x=round(tx,1),y=ty,`text-anchor`="middle",
              style="font-family:monospace;font-size:10px;fill:#4a8c68;",as.character(v)))}))
        )
      }

      div(class="viz-panel",style="padding:18px 14px;gap:18px;",
        div(
          div(style="display:flex;justify-content:space-between;align-items:baseline;margin-bottom:8px;",
            div(style="font-size:11px;font-weight:700;letter-spacing:.10em;text-transform:uppercase;color:#ff8096;",
                paste0("Field A  \u2014  ",nA," sample",if(nA==1)""else"s")),
            div(style="font-size:10px;color:var(--tx3);font-family:var(--mono);","90% interval")),
          make_iv_svg(a1,b1,"#ff8096","A")
        ),
        div(style="height:1px;background:var(--bdr);"),
        div(
          div(style="display:flex;justify-content:space-between;align-items:baseline;margin-bottom:8px;",
            div(style="font-size:11px;font-weight:700;letter-spacing:.10em;text-transform:uppercase;color:#5bc4f0;",
                paste0("Field B  \u2014  ",nB," sample",if(nB==1)""else"s")),
            div(style="font-size:10px;color:var(--tx3);font-family:var(--mono);","90% interval")),
          make_iv_svg(a2,b2,"#5bc4f0","B")
        )
      )

    } else if (mode == "gauge") {
      a1<-rv$a1;b1<-rv$b1;a2<-rv$a2;b2<-rv$b2;nA<-rv$nA;nB<-rv$nB
      mnA<-a1/(a1+b1);sdA<-sqrt(posterior_var(a1,b1))
      mnB<-a2/(a2+b2);sdB<-sqrt(posterior_var(a2,b2))

      make_gauge<-function(mn,sd,label,cbr,carc){
        # Layout (all px in viewBox units):
        # viewBox: 340 x 230
        # Field label: y=20 (top centre)
        # Semicircle: cx=170, cy=172, r=124
        #   arc top  = cy-r = 48  → clear of label at y=20
        #   arc base = cy   = 172
        # Tick labels sit at r+28 from centre, clamped so bottom labels are at most y=cy+10=182
        # Mean value: y=cy-28=144 (well inside the dial, above hub)
        # Sample count: y=cy+22=194 (below arc base, no overlap)
        # Total SVG height 230 gives 36px below cy for sample count
        cx<-170;cy<-172;r<-124;sw<-340;sh<-230
        p2d<-function(p)180-p*180
        arc_d<-function(d1,d2,rad){
          x1<-cx+rad*cos(d1*pi/180);y1<-cy-rad*sin(d1*pi/180)
          x2<-cx+rad*cos(d2*pi/180);y2<-cy-rad*sin(d2*pi/180)
          lf<-if(abs(d2-d1)>180)1 else 0
          paste0("M",round(x1,1),",",round(y1,1)," A",rad,",",rad," 0 ",lf,",0 ",round(x2,1),",",round(y2,1))
        }
        lo<-max(0,mn-sd);hi<-min(1,mn+sd)
        nd<-p2d(mn);nx<-cx+r*cos(nd*pi/180);ny<-cy-r*sin(nd*pi/180)

        # Major ticks at 0,0.1,...,1.0 — label positioned outside arc
        # For angles near 0° or 180° (far left/right), labels go outward horizontally
        # For angles near 90° (top), labels go above
        maj<-tagList(lapply(seq(0,1,.1),function(v){
          td<-p2d(v)
          rad_inner<-r-20;rad_outer<-r+10;rad_lbl<-r+26
          x1t<-cx+rad_inner*cos(td*pi/180);y1t<-cy-rad_inner*sin(td*pi/180)
          x2t<-cx+rad_outer*cos(td*pi/180);y2t<-cy-rad_outer*sin(td*pi/180)
          lx <-cx+rad_lbl  *cos(td*pi/180);ly <-cy-rad_lbl  *sin(td*pi/180)
          # Clamp label y so it doesn't go below the SVG
          ly_c<-min(ly+4, sh-8)
          tagList(
            tags$line(x1=round(x1t,1),y1=round(y1t,1),x2=round(x2t,1),y2=round(y2t,1),
              stroke="#5e9474",`stroke-width`="2",`stroke-linecap`="round"),
            tags$text(x=round(lx,1),y=round(ly_c,1),`text-anchor`="middle",
              `dominant-baseline`="middle",
              style="font-family:monospace;font-size:9.5px;fill:#5e9474;",as.character(v))
          )
        }))

        # Minor ticks at 0.05,0.15,...,0.95
        min_t<-tagList(lapply(seq(.05,.95,.1),function(v){
          td<-p2d(v)
          x1t<-cx+(r-12)*cos(td*pi/180);y1t<-cy-(r-12)*sin(td*pi/180)
          x2t<-cx+(r+4) *cos(td*pi/180);y2t<-cy-(r+4) *sin(td*pi/180)
          tags$line(x1=round(x1t,1),y1=round(y1t,1),x2=round(x2t,1),y2=round(y2t,1),
            stroke="#22402e",`stroke-width`="1.2",`stroke-linecap`="round")
        }))

        # Gradient defs for arc glow effect
        defs<-tags$defs(
          tags$radialGradient(id=paste0("hg_",gsub(" ","",label)),
            cx="50%",cy="100%",r="50%",
            tags$stop(offset="0%",  style=paste0("stop-color:",carc,";stop-opacity:0.18")),
            tags$stop(offset="100%",style="stop-color:#060f0a;stop-opacity:0")
          )
        )

        tags$svg(xmlns="http://www.w3.org/2000/svg",viewBox=paste0("0 0 ",sw," ",sh),
          style="width:100%;display:block;",
          defs,
          # Subtle glow fill under the arc area
          tags$path(d=paste0("M",cx-r-30,",",cy," Q",cx,",",cy-r-20," ",cx+r+30,",",cy," Z"),
            fill=paste0("url(#hg_",gsub(" ","",label),")"),opacity="0.6"),
          # Track arc
          tags$path(d=arc_d(180,0,r),fill="none",stroke="#0d1a10",`stroke-width`="24",`stroke-linecap`="round"),
          # Uncertainty arc
          tags$path(d=arc_d(p2d(hi),p2d(lo),r),fill="none",stroke=carc,`stroke-width`="24",`stroke-linecap`="round",opacity="0.85"),
          # Ticks
          min_t,maj,
          # Needle shadow
          tags$line(x1=cx+1.5,y1=cy+1.5,
            x2=round(cx+(nx-cx)*.9,1),y2=round(cy+(ny-cy)*.9,1),
            stroke="rgba(0,0,0,.4)",`stroke-width`="5",`stroke-linecap`="round"),
          # Needle
          tags$line(x1=cx,y1=cy,x2=round(nx,1),y2=round(ny,1),
            stroke="#ffffff",`stroke-width`="3",`stroke-linecap`="round"),
          # Hub
          tags$circle(cx=cx,cy=cy,r="12",fill=cbr),
          tags$circle(cx=cx,cy=cy,r="6", fill="#060f0a"),
          # Mean value — inside dial, well above hub and clear of arc bottom
          tags$text(x=cx,y=cy-32,`text-anchor`="middle",
            style=paste0("font-family:monospace;font-size:18px;fill:",cbr,";font-weight:700;"),
            sprintf("%.3f",mn)),
          # Field label — top of SVG, clear of arc top at y=48
          tags$text(x=cx,y=20,`text-anchor`="middle",
            style=paste0("font-family:monospace;font-size:12px;fill:",cbr,
                         ";font-weight:700;letter-spacing:.10em;"),label),
          # Sample count — below arc base, won't overlap tick labels
          tags$text(x=cx,y=cy+22,`text-anchor`="middle",
            style="font-family:monospace;font-size:11px;fill:#4a8c68;",
            paste0(if(label=="FIELD A")nA else nB," sample",
                   if((if(label=="FIELD A")nA else nB)==1)""else"s"))
        )
      }

      div(class="viz-panel",style="justify-content:center;align-items:center;padding:8px 0;",
        div(style="display:flex;gap:4px;justify-content:space-around;align-items:flex-start;width:100%;",
          make_gauge(mnA,sdA,"FIELD A","#ff8096","#e8404e"),
          make_gauge(mnB,sdB,"FIELD B","#5bc4f0","#2878c8")
        )
      )

    } else if (mode == "histogram") {
      a1<-rv$a1;b1<-rv$b1;a2<-rv$a2;b2<-rv$b2;nA<-rv$nA;nB<-rv$nB

      if (rv$hist_type == "silo") {
        # ── Per-silo probability histogram ─────────────────────
        ns <- rv$n_silos
        mnA <- a1/(a1+b1); mnB <- a2/(a2+b2)
        loA <- qbeta(0.05,a1,b1); hiA <- qbeta(0.95,a1,b1)
        loB <- qbeta(0.05,a2,b2); hiB <- qbeta(0.95,a2,b2)

        make_silo_hist <- function(mn, lo, hi, col, label, n_samp, uid) {
          sw<-340; sh<-160; pl<-28; pr<-16; pt<-28; pb<-30
          pw<-sw-pl-pr; ph<-sh-pt-pb; bw <- pw/ns; ay <- pt+ph
          grad_id <- paste0("shg",uid)
          defs <- tags$defs(tags$linearGradient(id=grad_id,x1="0",y1="0",x2="0",y2="1",
            tags$stop(offset="0%",  style=paste0("stop-color:",col,";stop-opacity:0.90")),
            tags$stop(offset="100%",style=paste0("stop-color:",col,";stop-opacity:0.22"))))
          grid_els <- tagList(lapply(c(0.25,0.5,0.75,1.0),function(v){
            gy<-pt+ph*(1-v)
            tags$line(x1=pl,y1=round(gy,1),x2=sw-pr,y2=round(gy,1),
              stroke=if(v==0.5)"#2a5040"else"#0f1e16",
              `stroke-width`=if(v==0.5)"1.2"else"0.6",
              `stroke-dasharray`=if(v==0.5)"4,3"else"none")}))
          band_top<-pt+ph*(1-hi); band_bot<-pt+ph*(1-lo)
          band<-tags$rect(x=pl,y=round(band_top,1),width=pw,
                          height=round(band_bot-band_top,1),fill=col,opacity="0.08",rx="2")
          bar_h<-ph*mn; bar_y<-pt+ph-bar_h
          bars<-tagList(lapply(seq_len(ns),function(i){
            bx<-pl+(i-1)*bw
            tagList(
              tags$rect(x=round(bx+1,1),y=round(bar_y,1),width=round(bw-2,1),height=round(bar_h,1),
                        rx="2",fill=paste0("url(#",grad_id,")"),opacity="0.88"),
              tags$rect(x=round(bx+1,1),y=round(bar_y,1),width=round(bw-2,1),height=round(bar_h,1),
                        rx="2",fill="none",stroke=col,`stroke-width`="0.5",opacity="0.4"))}))
          mean_y<-pt+ph*(1-mn)
          mean_line<-tags$line(x1=pl,y1=round(mean_y,1),x2=sw-pr,y2=round(mean_y,1),
                               stroke=col,`stroke-width`="1.5",`stroke-dasharray`="5,3",opacity="0.9")
          mean_lbl<-tags$text(x=sw-pr-2,y=round(mean_y-3,1),`text-anchor`="end",
                              style=paste0("font-family:monospace;font-size:10px;fill:",col,";font-weight:700;"),
                              sprintf("%.3f",mn))
          y_labels<-tagList(lapply(c(0,0.5,1.0),function(v){
            gy<-pt+ph*(1-v)
            tagList(tags$line(x1=pl-3,y1=round(gy,1),x2=pl,y2=round(gy,1),stroke="#22402e",`stroke-width`="1"),
                    tags$text(x=pl-5,y=round(gy+4,1),`text-anchor`="end",
                      style="font-family:monospace;font-size:9px;fill:#4a8c68;",as.character(v)))}))
          step<-if(ns<=10)1 else if(ns<=20)2 else 5
          x_labels<-tagList(lapply(seq(1,ns,by=step),function(i){
            tx<-pl+(i-0.5)*bw
            tags$text(x=round(tx,1),y=ay+14,`text-anchor`="middle",
                      style="font-family:monospace;font-size:9px;fill:#4a8c68;",as.character(i))}))
          tags$svg(xmlns="http://www.w3.org/2000/svg",viewBox=paste0("0 0 ",sw," ",sh),style="width:100%;display:block;",
            defs,tags$rect(x=pl,y=pt,width=pw,height=ph,fill="#030809",rx="3",opacity="0.7"),
            grid_els,band,bars,mean_line,mean_lbl,
            tags$line(x1=pl,y1=ay,x2=sw-pr,y2=ay,stroke="#22402e",`stroke-width`="1.2"),
            tags$line(x1=pl,y1=pt,x2=pl,y2=ay,stroke="#22402e",`stroke-width`="1"),
            y_labels,x_labels,
            tags$text(x=pl+2,y=18,`text-anchor`="start",
              style=paste0("font-family:monospace;font-size:11px;fill:",col,";font-weight:700;letter-spacing:.08em;"),label),
            tags$text(x=sw-pr-2,y=18,`text-anchor`="end",
              style="font-family:monospace;font-size:11px;fill:#4a8c68;",
              paste0(n_samp," sample",if(n_samp==1)""else"s")),
            tags$text(x=pl+pw/2,y=sh-2,`text-anchor`="middle",
              style="font-family:monospace;font-size:9px;fill:#3d6b50;","Silo number"))
        }
        div(class="viz-panel",style="padding:14px 10px;gap:14px;",
          make_silo_hist(mnA,loA,hiA,"#ff8096","FIELD A",nA,"A"),
          div(style="height:1px;background:var(--bdr);"),
          make_silo_hist(mnB,loB,hiB,"#5bc4f0","FIELD B",nB,"B"))

      } else {
        # density histogram — shows the beta posterior as bars
        nb<-rv$hist_bins
        mids<-seq(1/(2*nb),1-1/(2*nb),1/nb)
        dA<-dbeta(mids,a1,b1); dB<-dbeta(mids,a2,b2)
        dA[!is.finite(dA)]<-0; dB[!is.finite(dB)]<-0

        make_dens_hist<-function(dens,cf,label,ns,uid){
          sw<-340;sh<-106;pl<-22;pr<-14;pt<-22;pb<-24
          pw<-sw-pl-pr;ph<-sh-pt-pb;mx<-max(dens,.01);bw<-pw/nb;ay<-pt+ph
          grad_id<-paste0("hg",uid)
          defs<-tags$defs(tags$linearGradient(id=grad_id,x1="0",y1="0",x2="0",y2="1",
            tags$stop(offset="0%",  style=paste0("stop-color:",cf,";stop-opacity:0.92")),
            tags$stop(offset="100%",style=paste0("stop-color:",cf,";stop-opacity:0.28"))))
          grid_els<-tagList(lapply(c(.33,.66,1.0),function(f){
            gy<-pt+ph*(1-f)
            tags$line(x1=pl,y1=round(gy,1),x2=sw-pr,y2=round(gy,1),
              stroke=if(f==1.0)"#22402e"else"#0f1e16",
              `stroke-width`=if(f==1.0)"1.2"else"0.7")}))
          bars<-tagList(lapply(seq_len(nb),function(i){
            bx<-pl+(i-1)*bw;bh<-max(1,(dens[i]/mx)*ph);by<-pt+ph-bh
            tagList(
              tags$rect(x=round(bx+.8,1),y=round(by,1),width=round(bw-1.6,1),height=round(bh,1),
                        rx="2",fill=paste0("url(#",grad_id,")"),opacity="0.88"),
              tags$rect(x=round(bx+.8,1),y=round(by,1),width=round(bw-1.6,1),height=round(bh,1),
                        rx="2",fill="none",stroke=cf,`stroke-width`="0.6",opacity="0.35"))}))
          tags$svg(xmlns="http://www.w3.org/2000/svg",viewBox=paste0("0 0 ",sw," ",sh),style="width:100%;display:block;",
            defs,tags$rect(x=pl,y=pt,width=pw,height=ph,fill="#030809",rx="3",opacity="0.7"),
            grid_els,bars,
            tags$line(x1=pl,y1=ay,x2=sw-pr,y2=ay,stroke="#22402e",`stroke-width`="1.2"),
            tagList(lapply(c(0,.25,.5,.75,1),function(v){tx<-pl+v*pw;tagList(
              tags$line(x1=round(tx,1),y1=ay,x2=round(tx,1),y2=ay+3,stroke="#22402e",`stroke-width`="1"),
              tags$text(x=round(tx,1),y=ay+16,`text-anchor`="middle",
                style="font-family:monospace;font-size:9px;fill:#4a8c68;",as.character(v)))})),
            tags$text(x=pl+2,y=16,`text-anchor`="start",
              style=paste0("font-family:monospace;font-size:10px;fill:",cf,";font-weight:700;letter-spacing:.08em;"),label),
            tags$text(x=sw-pr-2,y=16,`text-anchor`="end",
              style="font-family:monospace;font-size:10px;fill:#4a8c68;",
              paste0(ns," sample",if(ns==1)""else"s")))
        }
        div(class="viz-panel",style="padding:14px 10px;gap:14px;",
          make_dens_hist(dA,"#ff8096","FIELD A",nA,"A"),
          div(style="height:1px;background:var(--bdr);"),
          make_dens_hist(dB,"#5bc4f0","FIELD B",nB,"B"))
      }
    }
  })

  output$plot <- renderPlot({
    req(rv$viz_mode=="curve")
    op<-par(no.readonly=TRUE);on.exit(par(op),add=TRUE)
    bg_col <- "#060f0a"
    par(mar=c(4.0,4.5,0.8,1.2),bg=bg_col,
        col.axis="#4a8c68",col.lab="#7ec49a",
        cex.axis=0.95,cex.lab=1.05,fg="#1a3024")
    x  <- seq(.001,.999,length.out=1000)
    y1 <- dbeta(x,rv$a1,rv$b1); y2 <- dbeta(x,rv$a2,rv$b2)
    ym <- max(y1,y2,2)*1.18
    plot.new(); plot.window(xlim=c(0,1),ylim=c(0,ym))
    # Darker background rect
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],
         col="#050d08",border=NA)
    # Subtle grid
    abline(h=pretty(c(0,ym),n=5),col="#0d1a10",lwd=0.8)
    abline(v=seq(0,1,.2),col="#0d1a10",lwd=0.7)
    abline(h=1,col=adjustcolor("#22402e",alpha.f=.5),lty=2,lwd=0.9)
    # Clean fill under each curve — single pass, natural shape
    polygon(c(x,rev(x)), c(y1,rep(0,length(x))),
            col=adjustcolor("#ff8096", alpha.f=.11), border=NA)
    polygon(c(x,rev(x)), c(y2,rep(0,length(x))),
            col=adjustcolor("#5bc4f0", alpha.f=.11), border=NA)
    # Soft outer glow then crisp line on top
    lines(x, y1, col=adjustcolor("#ff8096", alpha.f=.20), lwd=9)
    lines(x, y2, col=adjustcolor("#5bc4f0", alpha.f=.20), lwd=9)
    lines(x, y1, col="#ff8096", lwd=2.2)
    lines(x, y2, col="#5bc4f0", lwd=2.2)
    # Labels — 11px equivalent, clear of curves
    pAm <- rv$a1/(rv$a1+rv$b1); pBm <- rv$a2/(rv$a2+rv$b2)
    yA  <- min(dbeta(pAm,rv$a1,rv$b1)*.40,ym*.70)
    yB  <- min(dbeta(pBm,rv$a2,rv$b2)*.40,ym*.70)
    text(pAm,yA,sprintf("A\n%.3f",pAm),col="#ff8096",cex=0.88,font=2)
    text(pBm,yB,sprintf("B\n%.3f",pBm),col="#5bc4f0",cex=0.88,font=2)
    axis(1,col=NA,col.ticks="#1a3024",col.axis="#4a8c68")
    axis(2,col=NA,col.ticks="#1a3024",col.axis="#4a8c68",las=1)
    title(xlab="p  (missile probability per silo)",
          ylab="Posterior density",col.lab="#7ec49a")
    legend("topright",
           legend=c("Field A","Field B","Uniform prior"),
           col=c("#ff8096","#5bc4f0",adjustcolor("#d8efe2",alpha.f=.18)),
           lty=c(1,1,2),lwd=c(2.4,2.4,1.1),bty="n",
           text.col="#7ec49a",cex=0.84,bg=NA,inset=c(.01,.01))
  },bg="#060f0a")

  # rules popup
  observeEvent(input$rules,{
    showModal(modalDialog(title="How It Works",easyClose=TRUE,footer=modalButton("Close"),
      HTML(paste0(
        "<div style='font-size:14px;line-height:1.9;color:#9ecfb5;'>",
        "<p><strong style='color:#d8efe2'>Goal:</strong> Two fields each contain a hidden number of missiles. ",
        "You have ",rv$start_units," coverage units to split between them. Allocate enough to each field to block all its missiles.</p>",
        "<p><strong style='color:#d8efe2'>Sampling:</strong> Press <em>Sample Both Fields</em> to take one reading from each field. ",
        "This updates the probability estimates. However, every combined sample also rolls against a rising detonation risk.</p>",
        "<p><strong style='color:#d8efe2'>Detonation risk:</strong> Starts at 0%. Rises with every combined sample. ",
        "If the roll fires, the trial ends immediately and you lose.</p>",
        "<p><strong style='color:#d8efe2'>Practice trial:</strong> The first round is practice. Results are not recorded.</p>",
        "<p><strong style='color:#d8efe2'>After ",N_TRIALS," real trials</strong> the study ends and your results are shown.</p>",
        "</div>"
      ))
    ))
  })

  # tutorial screen — shown before practice when enabled
  output$tutorialScreen <- renderUI({
    if (rv$screen != "tutorial") return(NULL)
    idx <- rv$tutorial_index
    if (idx < 1 || idx > length(TUTORIAL_SCENARIOS)) return(NULL)
    sc  <- TUTORIAL_SCENARIOS[[idx]]
    U   <- rv$start_units

    div(class="consent-screen",
      div(class="consent-card", style="max-width:680px;",
        # Header
        div(style="display:flex;align-items:center;gap:12px;margin-bottom:6px;",
            div(
            div(style="font-size:11px;font-weight:700;letter-spacing:.14em;text-transform:uppercase;color:var(--tx3);",
                sc$label),
            div(style="font-size:18px;font-weight:800;color:var(--tx);","Allocation Practice")
          )
        ),
        div(class="consent-divider"),

        # Explanation
        div(style="font-size:13px;color:#9ecfb5;line-height:1.8;margin-bottom:18px;",
          "In this round the field probabilities are already known. ",
          strong(style="color:var(--tx);","No sampling needed — just allocate your units and launch."),
          " Use the sliders below to split your ", strong(U), " units between the two fields."
        ),

        # Hint box
        div(style="background:var(--accb);border:1px solid var(--accm);border-radius:var(--r2);padding:10px 14px;margin-bottom:20px;font-size:12px;color:#7ec49a;",
          sc$hint
        ),

        # Probability display
        div(style="display:grid;grid-template-columns:1fr 1fr;gap:10px;margin-bottom:20px;",
          div(class="field-row red", style="flex-direction:column;align-items:flex-start;gap:6px;",
            span(class="fr-tag red","Field A"),
            div(style="font-family:var(--mono);font-size:28px;font-weight:800;color:#ff8096;",
                sprintf("%.0f%%", sc$pA * 100)),
            div(style="font-size:11px;color:var(--tx2);",
                paste0("per-silo missile probability"))
          ),
          div(class="field-row blue", style="flex-direction:column;align-items:flex-start;gap:6px;",
            span(class="fr-tag blue","Field B"),
            div(style="font-family:var(--mono);font-size:28px;font-weight:800;color:#5bc4f0;",
                sprintf("%.0f%%", sc$pB * 100)),
            div(style="font-size:11px;color:var(--tx2);",
                paste0("per-silo missile probability"))
          )
        ),

        # Allocation sliders
        div(style="margin-bottom:6px;",
          div(style="display:flex;justify-content:space-between;align-items:baseline;margin-bottom:4px;",
            span(style="font-size:10px;font-weight:700;letter-spacing:.10em;text-transform:uppercase;color:#ff8096;","Field A"),
            span(style="font-family:var(--mono);font-size:16px;font-weight:700;color:#ff8096;",
                 textOutput("tut_uA_display", inline=TRUE))
          ),
          sliderInput("tut_uA", NULL, min=0, max=U, value=floor(U/2), step=1, width="100%")
        ),
        div(style="margin-bottom:20px;",
          div(style="display:flex;justify-content:space-between;align-items:baseline;margin-bottom:4px;",
            span(style="font-size:10px;font-weight:700;letter-spacing:.10em;text-transform:uppercase;color:#5bc4f0;","Field B"),
            span(style="font-family:var(--mono);font-size:16px;font-weight:700;color:#5bc4f0;",
                 textOutput("tut_uB_display", inline=TRUE))
          ),
          sliderInput("tut_uB", NULL, min=0, max=U, value=ceiling(U/2), step=1, width="100%")
        ),

        actionButton("tut_launch", "Launch \u2192", class="btn-primary", width="100%")
      )
    )
  })

  # keep tutorial sliders linked
  observeEvent(input$tut_uA, {
    req(!is.null(input$tut_uA))
    new_uB <- rv$start_units - as.integer(input$tut_uA)
    if (!is.null(input$tut_uB) && as.integer(input$tut_uB) != new_uB)
      updateSliderInput(session, "tut_uB", value=new_uB)
  }, ignoreInit=TRUE)

  observeEvent(input$tut_uB, {
    req(!is.null(input$tut_uB))
    new_uA <- rv$start_units - as.integer(input$tut_uB)
    if (!is.null(input$tut_uA) && as.integer(input$tut_uA) != new_uA)
      updateSliderInput(session, "tut_uA", value=new_uA)
  }, ignoreInit=TRUE)

  output$tut_uA_display <- renderText(as.character(input$tut_uA %||% floor(rv$start_units/2)))
  output$tut_uB_display <- renderText(as.character(rv$start_units - (input$tut_uA %||% floor(rv$start_units/2))))

  # score the tutorial round and show result
  observeEvent(input$tut_launch, {
    req(rv$screen == "tutorial")
    idx <- rv$tutorial_index
    sc  <- TUTORIAL_SCENARIOS[[idx]]
    U   <- rv$start_units
    uA  <- min(as.integer(input$tut_uA %||% floor(U/2)), U)
    uB  <- U - uA

    hitA  <- min(sc$mA, uA); hitB <- min(sc$mB, uB)
    missA <- max(sc$mA - uA, 0); missB <- max(sc$mB - uB, 0)
    won   <- (missA + missB) == 0

    # get bayesian optimal for display
    opt <- bayes_optimal_split(
      sc$pA * 10 + 1, (1 - sc$pA) * 10 + 1,
      sc$pB * 10 + 1, (1 - sc$pB) * 10 + 1,
      rv$n_silos, U)

    icon_s  <- if(won) '<span style="font-size:11px;font-weight:800;letter-spacing:.14em;color:#3db87a;">ALL COVERED</span>' else '<span style="font-size:11px;font-weight:800;letter-spacing:.14em;color:#c96070;">MISSILES ESCAPED</span>'
    col_s   <- if(won) "#3db87a" else "#c96070"
    result_txt <- if(won) "All missiles covered!" else
      paste0(missA + missB, " missile(s) escaped.")
    split_txt <- if(!is.na(opt$opt_uA))
      sprintf("You: A%d/B%d &nbsp;&middot;&nbsp; Optimal: A%d/B%d", uA, uB, opt$opt_uA, opt$opt_uB)
      else sprintf("You: A%d/B%d", uA, uB)

    is_last <- idx == length(TUTORIAL_SCENARIOS)
    next_label <- if(is_last) "Continue to Practice Round \u2192" else paste0("Next Tutorial \u2192 (", idx+1, "/", length(TUTORIAL_SCENARIOS), ")")

    showModal(modalDialog(
      easyClose = FALSE,
      footer = actionButton("tut_next", next_label, class="btn-primary"),
      tags$div(
        style=paste0("border-top:3px solid ",col_s,";margin:-1px -1px 0;border-radius:var(--r3) var(--r3) 0 0;"),
        div(style="text-align:center;padding:16px 0 8px;",
          div(style="margin-bottom:12px;", HTML(icon_s)),
          div(style="font-size:16px;font-weight:700;color:var(--tx);margin-bottom:6px;", result_txt),
          div(style="font-size:12px;color:var(--tx2);font-family:var(--mono);margin-bottom:6px;",
              HTML(split_txt)),
          div(style="font-size:11px;color:var(--tx3);", sc$label)
        )
      )
    ))
  })

  # move to next tutorial or into game
  observeEvent(input$tut_next, {
    removeModal()
    if (rv$tutorial_index < length(TUTORIAL_SCENARIOS)) {
      rv$tutorial_index <- rv$tutorial_index + 1
      # Reset sliders for next tutorial
      U <- rv$start_units
      updateSliderInput(session, "tut_uA", min=0, max=U, value=floor(U/2))
      updateSliderInput(session, "tut_uB", min=0, max=U, value=ceiling(U/2))
    } else {
      rv$tutorial_done <- TRUE
      rv$screen <- "game"
      start_trial()
    }
  })

  # results summary shown after all trials
  output$resultsScreen <- renderUI({
    if (rv$screen!="results") return(NULL)
    h<-rv$trial_history;n<-length(h)
    wins    <- sum(sapply(h,function(x)isTRUE(x$win)))
    dets    <- sum(sapply(h,function(x)isTRUE(x$detected)))
    win_pct <- round(wins/max(n,1)*100)
    avg_samp<- if(n>0)round(mean(sapply(h,function(x)x$samples)),1) else 0
    avg_time<- if(n>0)round(mean(sapply(h,function(x)x$elapsed)),1) else 0
    errs    <- Filter(function(x)!is.na(x),sapply(h,function(x)x$split_err))
    avg_err <- if(length(errs)>0)round(mean(errs),2) else NA_real_

    vl <- switch(rv$viz_mode,"off"="None","curve"="Curve","dotplot"="Interval Bar",
                 "gauge"="Gauge","histogram"="Histogram","None")

    make_row<-function(hh,delay){
      icon  <-if(isTRUE(hh$win))"WIN"else if(isTRUE(hh$detected))"DET"else"FAIL"
      cls   <-if(isTRUE(hh$win))"rs-trial-row win"else if(isTRUE(hh$detected))"rs-trial-row det"else"rs-trial-row lose"
      label <-if(isTRUE(hh$win))"Success"else if(isTRUE(hh$detected))"Detonated"else"Failed"
      # True optimal = units needed to exactly cover each field's actual missiles
      err_t <- if(!is.na(hh$split_err) && !is.null(hh$true_opt_A) && !is.na(hh$true_opt_A)){
        paste0("A:",hh$player_uA," true:A",hh$true_opt_A,
               " (",if(hh$split_err==0)"on target"else paste0(hh$split_err,"u off"),")")
      } else "\u2014"
      div(class=cls,style=paste0("animation-delay:",delay,"ms;"),
        span(class="rs-trial-icon",icon),span(class="rs-trial-num",paste0("Trial ",hh$trial)),
        span(class="rs-trial-label",label),span(class="rs-trial-stat",paste0(hh$samples," samples")),
        span(class="rs-trial-stat",paste0(hh$elapsed,"s")),span(class="rs-trial-stat err",err_t))
    }

    div(class="results-screen",
      div(class="rs-header",
        div(class="rs-title","Study Complete"),
        div(class="rs-subtitle",paste0("Code: ",rv$game_code,"  \u00b7  ",N_TRIALS," trials  \u00b7  Risk Mode"))
      ),
      div(class="rs-section-title","Session Parameters"),
      div(class="rs-settings-grid",
        div(class="rs-setting",div(class="rs-set-k","Visualization"),div(class="rs-set-v",vl)),
        div(class="rs-setting",div(class="rs-set-k","Coverage Units"),div(class="rs-set-v",rv$start_units)),
        div(class="rs-setting",div(class="rs-set-k","Total Missiles"),div(class="rs-set-v",rv$total_missiles)),
        div(class="rs-setting",div(class="rs-set-k","Data Collected"),div(class="rs-set-v",if(rv$consented)"Yes"else"No (declined)")),
        div(class="rs-setting",div(class="rs-set-k","Risk k"),div(class="rs-set-v",rv$risk_k)),
        div(class="rs-setting",div(class="rs-set-k","Risk exp"),div(class="rs-set-v",rv$risk_exp))
      ),
      div(class="rs-section-title","Performance Summary"),
      div(class="rs-stat-grid",
        div(class="rs-stat-card",style="animation-delay:50ms;",div(class="rs-stat-val green",paste0(wins,"/",n)),div(class="rs-stat-label","Missions\nSucceeded")),
        div(class="rs-stat-card",style="animation-delay:120ms;",div(class="rs-stat-val",paste0(win_pct,"%")),div(class="rs-stat-label","Success\nRate")),
        div(class="rs-stat-card",style="animation-delay:190ms;",div(class="rs-stat-val amber",paste0(dets,"/",n)),div(class="rs-stat-label","Detonations")),
        div(class="rs-stat-card",style="animation-delay:260ms;",div(class="rs-stat-val",avg_samp),div(class="rs-stat-label","Avg Samples\nPer Trial")),
        div(class="rs-stat-card",style="animation-delay:330ms;",div(class="rs-stat-val",if(!is.na(avg_err))avg_err else"\u2014"),div(class="rs-stat-label","Avg Allocation\nError"))
      ),
      div(class="rs-section-title","Trial Breakdown"),
      div(class="rs-trial-list",
          tagList(lapply(seq_along(h),function(i)make_row(h[[i]],100+i*60)))),
      if(rv$is_admin){tagList(
        div(class="rs-section-title","Export Data"),
        div(class="rs-dl-row",
          div(class="rs-dl-card",
            div(class="rs-dl-label","All Events"),
            div(class="rs-dl-desc","Every sample and launch event logged"),
            downloadButton("dl_events","Download CSV",class="btn-primary rs-dl-btn")),
          div(class="rs-dl-card",
            div(class="rs-dl-label","Launch Trials"),
            div(class="rs-dl-desc","One row per trial with all study variables"),
            downloadButton("dl_launches","Download CSV",class="btn-primary rs-dl-btn"))
        )
      )}
    )
  })

  output$dl_events <- downloadHandler(
    filename=function()paste0("events_",Sys.Date(),".csv"),
    content=function(f){if(file.exists(EVENTS_FILE))file.copy(EVENTS_FILE,f)else write.csv(data.frame(msg="No data"),f,row.names=FALSE)}
  )
  output$dl_launches <- downloadHandler(
    filename=function()paste0("launches_",Sys.Date(),".csv"),
    content=function(f){if(file.exists(LAUNCHES_FILE))file.copy(LAUNCHES_FILE,f)else write.csv(data.frame(msg="No data"),f,row.names=FALSE)}
  )
}

shinyApp(ui, server)

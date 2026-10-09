// =============================================================================
//  include/yjk/post/ViewerTemplate.h  ——  后处理界面模板（单文件 HTML）
//
//  【为什么不用 three.js / WebGL】
//    · 结果文件要能离线双击打开，CDN 依赖在无网环境直接白屏
//    · 结构模型是【线框 + 少量四边形】，Canvas 2D 完全够（几千根杆 60fps）
//    · 不引入构建步骤：模板 + JSON 替换 → 一个文件
//
//  【三维投影】手写 yaw/pitch 旋转 + 画家算法深度排序。
//  不用 z-buffer 是因为结构模型的图元都是细长线段，
//  按中点深度排序在视觉上与真三维几乎不可分，而实现只要几行。
//
//  【为什么必须做深度排序】
//  不排序时后画的杆永远压在先画的杆上，旋转到背面会看到
//  "前面被后面遮挡"的错乱 —— 这是最容易被当成"求解错了"的界面 bug。
// =============================================================================
#pragma once

namespace yjk {
namespace post {

inline const char* kViewerTemplate = R"YJKRAW(<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__YJK_TITLE__</title>
<style>
  :root{
    --bg:#12161c; --panel:#1b2029; --panel2:#232a35; --line:#333c4a;
    --fg:#dfe4ec; --dim:#8b96a8; --accent:#4aa3ff; --ok:#37c76a; --warn:#ffb648;
    --bad:#ff5f56;
    --font:"Segoe UI","Microsoft YaHei",system-ui,sans-serif;
  }
  *{box-sizing:border-box}
  html,body{height:100%;margin:0}
  body{background:var(--bg);color:var(--fg);font-family:var(--font);font-size:13px;
       display:flex;flex-direction:column;overflow:hidden}
  header{display:flex;align-items:center;gap:14px;padding:8px 14px;background:var(--panel);
         border-bottom:1px solid var(--line);flex:0 0 auto;flex-wrap:wrap}
  header h1{font-size:14px;margin:0;font-weight:600;letter-spacing:.3px}
  header .sub{color:var(--dim);font-size:12px}
  .badge{padding:2px 8px;border-radius:10px;font-size:11px;background:var(--panel2);
         border:1px solid var(--line);color:var(--dim)}
  .badge.ok{color:var(--ok);border-color:#2b6b45}
  .badge.bad{color:var(--bad);border-color:#7a3230}
  main{flex:1 1 auto;display:flex;min-height:0}
  #side{width:258px;flex:0 0 258px;background:var(--panel);border-right:1px solid var(--line);
        overflow-y:auto;padding:10px}
  #stage{flex:1 1 auto;position:relative;min-width:0}
  canvas#gl{display:block;width:100%;height:100%;cursor:grab;background:#0d1116}
  canvas#gl.drag{cursor:grabbing}
  #right{width:300px;flex:0 0 300px;background:var(--panel);border-left:1px solid var(--line);
         overflow-y:auto;padding:10px;display:none}
  #right.on{display:block}
  .grp{margin-bottom:12px}
  .grp>h3{margin:0 0 6px;font-size:11px;color:var(--dim);text-transform:uppercase;
          letter-spacing:.8px;font-weight:600}
  label.row{display:flex;align-items:center;gap:6px;padding:3px 0;cursor:pointer;font-size:12px}
  select,input[type=range]{width:100%}
  select{background:var(--panel2);color:var(--fg);border:1px solid var(--line);
         border-radius:4px;padding:4px 6px;font-size:12px}
  button{background:var(--panel2);color:var(--fg);border:1px solid var(--line);
         border-radius:4px;padding:4px 9px;font-size:12px;cursor:pointer}
  button:hover{border-color:var(--accent);color:#fff}
  button.act{background:var(--accent);color:#05121f;border-color:var(--accent);font-weight:600}
  .btns{display:flex;gap:5px;flex-wrap:wrap}
  .kv{display:flex;justify-content:space-between;gap:8px;padding:2px 0;font-size:12px}
  .kv span:first-child{color:var(--dim)}
  .kv span:last-child{font-variant-numeric:tabular-nums}
  table{width:100%;border-collapse:collapse;font-size:11px;font-variant-numeric:tabular-nums}
  th,td{padding:3px 4px;text-align:right;border-bottom:1px solid var(--line)}
  th{color:var(--dim);font-weight:600;text-align:right}
  th:first-child,td:first-child{text-align:left}
  #legend{position:absolute;left:14px;bottom:14px;background:rgba(20,25,33,.92);
          border:1px solid var(--line);border-radius:6px;padding:8px 10px;pointer-events:none}
  #legend canvas{display:block;border-radius:3px}
  #legend .lab{display:flex;justify-content:space-between;font-size:10px;color:var(--dim);margin-top:3px}
  #hint{position:absolute;right:14px;top:12px;color:var(--dim);font-size:11px;
        background:rgba(20,25,33,.85);padding:5px 9px;border-radius:5px;pointer-events:none}
  #status{position:absolute;left:14px;top:12px;color:var(--dim);font-size:11px;
          background:rgba(20,25,33,.85);padding:5px 9px;border-radius:5px;pointer-events:none}
  .warnbox{background:#2a2018;border:1px solid #6b4a22;color:var(--warn);
           padding:6px 8px;border-radius:4px;font-size:11px;margin-bottom:10px}
  #pick{position:absolute;background:rgba(20,25,33,.96);border:1px solid var(--accent);
        border-radius:5px;padding:6px 9px;font-size:11px;pointer-events:none;display:none;
        z-index:5;white-space:pre}
</style>
</head>
<body>
<header>
  <h1>__YJK_TITLE__</h1>
  <span class="sub" id="hsub"></span>
  <span class="badge" id="hok"></span>
  <span class="badge" id="hres"></span>
  <span class="badge" id="hsize"></span>
</header>
<main>
  <div id="side">
    <div id="warnzone"></div>
    <div class="grp">
      <h3>视图</h3>
      <div class="btns">
        <button data-view="iso">等轴测</button>
        <button data-view="top">俯视</button>
        <button data-view="front">正视</button>
        <button data-view="side">侧视</button>
        <button id="fit">适应窗口</button>
      </div>
    </div>
    <div class="grp">
      <h3>着色</h3>
      <select id="cmode"></select>
      <div class="kv"><span>范围</span><span id="crange">—</span></div>
    </div>
    <div class="grp">
      <h3>变形</h3>
      <label class="row"><input type="checkbox" id="deform" checked>显示变形</label>
      <input type="range" id="dscale" min="0" max="1000" value="500">
      <div class="kv"><span>放大系数</span><span id="dscaletxt">1×</span></div>
      <div class="btns"><button id="play">▶ 动画</button><button id="reset">复位</button></div>
      <label class="row"><input type="checkbox" id="undef">显示原轮廓</label>
    </div>
    <div class="grp">
      <h3>显示</h3>
      <label class="row"><input type="checkbox" id="snode">节点</label>
      <label class="row"><input type="checkbox" id="sbeam" checked>梁 / 柱</label>
      <label class="row"><input type="checkbox" id="sshell" checked>板 / 壳</label>
      <label class="row"><input type="checkbox" id="sreact" checked>支座反力</label>
      <label class="row"><input type="checkbox" id="saxes" checked>坐标轴</label>
      <label class="row"><input type="checkbox" id="slabel">节点号</label>
    </div>
    <div class="grp">
      <h3>模型统计</h3>
      <div id="stats"></div>
    </div>
    <div class="grp">
      <h3>楼层指标</h3>
      <div id="stories"></div>
    </div>
  </div>
  <div id="stage">
    <canvas id="gl"></canvas>
    <div id="status"></div>
    <div id="hint">左键拖动旋转 · 滚轮缩放 · 右键平移 · 单击构件看内力</div>
    <div id="legend"><canvas id="lg" width="180" height="12"></canvas>
      <div class="lab"><span id="lgmin">—</span><span id="lgname">—</span><span id="lgmax">—</span></div>
    </div>
    <div id="pick"></div>
  </div>
  <div id="right">
    <div class="grp"><h3>构件信息</h3><div id="minfo"></div></div>
    <div class="grp"><h3>内力图（沿杆长）</h3>
      <canvas id="dg" width="276" height="320"></canvas>
      <div class="kv" style="margin-top:6px"><span>图例</span><span>实线=绕 z(Mz) / 虚线=绕 y(My)</span></div>
    </div>
    <div class="grp"><h3>沿杆极值</h3><div id="mstat"></div></div>
    <div class="btns"><button id="close">关闭</button></div>
  </div>
</main>
<script>
"use strict";
const DATA = __YJK_DATA__;

// ---- 数值工具 ----
const F = (v,d=4)=>{ if(!isFinite(v)) return "—";
  const a=Math.abs(v);
  if(a!==0 && (a<1e-3||a>=1e6)) return v.toExponential(2);
  return v.toFixed(d); };
const clamp=(v,a,b)=>v<a?a:(v>b?b:v);

// ---- 索引常量（与 C++ _schema 对应）----
const NX=0,NY=1,NZ=2,NSTORY=3,NMASK=4,NUX=5,NUY=6,NUZ=7,NRX=8,NRY=9,NRZ=10,NRX2=11,NRY2=12,NRZ2=13;
const SXI=0,SN=1,SVY=2,SVZ=3,ST=4,SMY=5,SMZ=6,SSMAX=7,SSMIN=8,SVM=9;
const BID=0,BNI=1,BNJ=2,BKIND=3,BL=4,BSMAX=5,BSMIN=6,BVM=7,BUTIL=8,BMMAX=9,BVMAX=10,BNMAX=11,BNMIN=12,BW=13,BST=14;
const QID=0,QN0=1,QN1=2,QN2=3,QN3=4,QAREA=5,QT=6,QSX=7,QSY=8,QTXY=9,QS1=10,QS2=11,QVM=12,QW=13;
const KINDNAME=DATA.memberKindNames||["柱","梁","支撑","其它"];

// ---- 状态 ----
const S={ yaw:-0.62, pitch:0.42, zoom:1, panX:0, panY:0,
          dispScale:1, deform:true, undef:false,
          cmode:"", sel:-1, anim:false, animT:0,
          snode:false, sbeam:true, sshell:true, sreact:true, sax:true, slabel:false };

const cv=document.getElementById("gl"), ctx=cv.getContext("2d");
const nodes=DATA.nodes||[], beams=DATA.beams||[], shells=DATA.shells||[];
const fields=DATA.fields||[], stories=DATA.stories||[], env=DATA.envelope||{}, meta=DATA.meta||{};

// ---- 包围盒 ----
const BB=(()=>{ let a=[1e30,1e30,1e30],b=[-1e30,-1e30,-1e30];
  for(const n of nodes){ a[0]=Math.min(a[0],n[NX]);b[0]=Math.max(b[0],n[NX]);
    a[1]=Math.min(a[1],n[NY]);b[1]=Math.max(b[1],n[NY]);
    a[2]=Math.min(a[2],n[NZ]);b[2]=Math.max(b[2],n[NZ]); }
  if(a[0]>b[0]){a=[0,0,0];b=[1,1,1];}
  return {min:a,max:b,cx:(a[0]+b[0])/2,cy:(a[1]+b[1])/2,cz:(a[2]+b[2])/2,
          size:Math.max(b[0]-a[0],b[1]-a[1],b[2]-a[2],1e-6)};})();

// ---- 颜色映射（蓝→青→绿→黄→红，经典有限元云图）----
function jet(t){
  t=clamp(t,0,1);
  const stops=[[0,0.23,0.60],[0,0.55,0.85],[0,0.80,0.72],[0.35,0.88,0.35],
               [0.98,0.92,0.20],[1,0.55,0.10],[0.92,0.16,0.12]];
  const x=t*(stops.length-1); const i=Math.min(Math.floor(x),stops.length-2); const f=x-i;
  const a=stops[i],b=stops[i+1];
  return [Math.round(255*(a[0]+(b[0]-a[0])*f)),
          Math.round(255*(a[1]+(b[1]-a[1])*f)),
          Math.round(255*(a[2]+(b[2]-a[2])*f))];
}
const rgb=c=>`rgb(${c[0]},${c[1]},${c[2]})`;

// ---- 当前着色场 ----
function curField(){
  if(S.cmode.startsWith("f:")){ const i=+S.cmode.slice(2); return fields[i]||null; }
  if(S.cmode==="vm")   return {name:"von Mises 应力",unit:"kPa",v:beams.map(b=>b[BVM]),per:"beam"};
  if(S.cmode==="util") return {name:"应力比",unit:"—",v:beams.map(b=>b[BUTIL]),per:"beam"};
  if(S.cmode==="axial")return {name:"轴力",unit:"kN",v:beams.map(b=>b[BNMAX]),per:"beam",signed:false};
  return null;
}
function fieldRange(f){
  if(!f) return [0,1];
  if(f.per==="beam"){
    let mn=1e30,mx=-1e30; for(const v of f.v){ if(v<mn)mn=v; if(v>mx)mx=v; }
    if(mn>mx) return [0,1];
    return [mn,mx];
  }
  return [f.min,f.max];
}

// ---- 投影 ----
let VP={cx:0,cy:0,s:1};
function project(x,y,z){
  const dx=x-BB.cx, dy=y-BB.cy, dz=z-BB.cz;
  const cy_=Math.cos(S.yaw), sy_=Math.sin(S.yaw);
  let X=dx*cy_-dy*sy_, Y=dx*sy_+dy*cy_, Z=dz;
  const cp=Math.cos(S.pitch), sp=Math.sin(S.pitch);
  let Y2=Y*cp-Z*sp, Z2=Y*sp+Z*cp;
  // 正交投影（结构模型用正交更符合工程习惯，尺寸可比）
  return [VP.cx+S.panX+X*VP.s, VP.cy+S.panY-Z2*VP.s, Y2];
}
function nodePos(i,withDisp){
  const n=nodes[i];
  let x=n[NX],y=n[NY],z=n[NZ];
  if(withDisp&&S.deform){ const k=S.dispScale*(S.anim?S.animT:1);
    x+=n[NUX]*k; y+=n[NUY]*k; z+=n[NUZ]*k; }
  return [x,y,z];
}

// ---- 绘制 ----
function resize(){
  const r=cv.parentElement.getBoundingClientRect();
  const dpr=window.devicePixelRatio||1;
  cv.width=Math.max(1,Math.round(r.width*dpr)); cv.height=Math.max(1,Math.round(r.height*dpr));
  ctx.setTransform(dpr,0,0,dpr,0,0);
  VP.cx=r.width/2; VP.cy=r.height/2;
  VP.s=(Math.min(r.width,r.height)/BB.size)*0.72*S.zoom;
  draw();
}
function draw(){
  const r=cv.getBoundingClientRect();
  ctx.clearRect(0,0,r.width,r.height);
  ctx.fillStyle="#0d1116"; ctx.fillRect(0,0,r.width,r.height);
  const f=curField(); const [fmin,fmax]=fieldRange(f);
  const fspan=(fmax-fmin)||1;

  // --- 壳（先画，作为背景）---
  if(S.sshell&&shells.length){
    const items=[];
    for(let qi=0;qi<shells.length;++qi){
      const q=shells[qi]; const P=[];
      let ok=true,dep=0;
      for(let k=0;k<4;++k){ const p=project(...nodePos(q[QN0+k],true));
        if(!isFinite(p[0])){ok=false;break;} P.push(p); dep+=p[2]; }
      if(!ok)continue; items.push({q,P,dep:dep/4});
    }
    items.sort((a,b)=>a.dep-b.dep);
    for(const it of items){
      const q=it.q;
      // 节点场时按四角平均；梁级场（von Mises / 应力比 / 轴力）对壳不适用，
      // 退化到壳自身的 von Mises，避免出现"着色与对象不匹配"的假云图。
      const v=(f&&f.per!=="beam") ? (f.v[q[QN0]]+f.v[q[QN1]]+f.v[q[QN2]]+f.v[q[QN3]])/4 : q[QVM];
      ctx.beginPath();
      ctx.moveTo(it.P[0][0],it.P[0][1]);
      for(let k=1;k<4;++k) ctx.lineTo(it.P[k][0],it.P[k][1]);
      ctx.closePath();
      ctx.fillStyle = f ? rgb(jet((v-fmin)/fspan)).replace("rgb","rgba").replace(")",",0.62)")
                        : "rgba(70,86,110,0.45)";
      ctx.fill();
      ctx.strokeStyle="rgba(150,170,200,0.30)"; ctx.lineWidth=0.6; ctx.stroke();
    }
  }

  // --- 梁 / 柱 ---
  if(S.sbeam){
    const items=[];
    for(let bi=0;bi<beams.length;++bi){
      const b=beams[bi];
      const a=project(...nodePos(b[BNI],true)), c=project(...nodePos(b[BNJ],true));
      if(!isFinite(a[0])||!isFinite(c[0]))continue;
      items.push({b,i:bi,a,c,dep:(a[2]+c[2])/2});
    }
    items.sort((x,y)=>x.dep-y.dep);
    for(const it of items){
      const b=it.b;
      const isSel=(S.sel===it.i);
      if(f&&f.per!=="beam"){
        // 节点场 → 沿杆渐变：分 6 段
        const va=f.v[b[BNI]], vb=f.v[b[BNJ]];
        const N=6;
        for(let k=0;k<N;++k){
          const t0=k/N,t1=(k+1)/N;
          const p0=project(...lerpNode(b[BNI],b[BNJ],t0)), p1=project(...lerpNode(b[BNI],b[BNJ],t1));
          const vm=(va+(vb-va)*((t0+t1)/2)-fmin)/fspan;
          ctx.beginPath(); ctx.moveTo(p0[0],p0[1]); ctx.lineTo(p1[0],p1[1]);
          ctx.strokeStyle=rgb(jet(vm));
          ctx.lineWidth=isSel?4.0:(b[BKIND]===0?3.0:2.2); ctx.stroke();
        }
      }else{
        const v=f?f.v[it.i]:0;
        ctx.beginPath(); ctx.moveTo(it.a[0],it.a[1]); ctx.lineTo(it.c[0],it.c[1]);
        ctx.strokeStyle=f?rgb(jet((v-fmin)/fspan)):(isSel?"#ffd45e":"#7f8ea6");
        ctx.lineWidth=isSel?4.0:(b[BKIND]===0?3.0:2.2);
        ctx.stroke();
      }
      if(isSel){ ctx.beginPath(); ctx.moveTo(it.a[0],it.a[1]); ctx.lineTo(it.c[0],it.c[1]);
        ctx.strokeStyle="rgba(255,212,94,0.35)"; ctx.lineWidth=9; ctx.stroke(); }
    }
    // 原轮廓
    if(S.undef&&S.deform){
      ctx.save(); ctx.setLineDash([3,3]); ctx.strokeStyle="#4c586b"; ctx.lineWidth=1;
      for(const b of beams){
        const a=project(...nodePos(b[BNI],false)), c=project(...nodePos(b[BNJ],false));
        ctx.beginPath(); ctx.moveTo(a[0],a[1]); ctx.lineTo(c[0],c[1]); ctx.stroke();
      }
      ctx.restore();
    }
  }

  // --- 支座反力箭头 ---
  if(S.sreact){
    let mx=0; for(const n of nodes) mx=Math.max(mx,Math.abs(n[NRX2]),Math.abs(n[NRY2]),Math.abs(n[NRZ2]));
    if(mx>1e-9){
      const L=BB.size*0.10;
      for(let i=0;i<nodes.length;++i){
        const n=nodes[i];
        const v=[n[NRX2],n[NRY2],n[NRZ2]];
        const mag=Math.hypot(v[0],v[1],v[2]);
        if(mag<1e-7*mx)continue;
        const p=project(...nodePos(i,true));
        const d=[v[0]/mx,v[1]/mx,v[2]/mx];
        const tip=project(...( (()=>{const q=nodePos(i,true);
          return [q[0]+d[0]*L,q[1]+d[1]*L,q[2]+d[2]*L];})() ));
        ctx.beginPath(); ctx.moveTo(p[0],p[1]); ctx.lineTo(tip[0],tip[1]);
        ctx.strokeStyle="#ff7a5e"; ctx.lineWidth=1.6; ctx.stroke();
        ctx.beginPath(); ctx.arc(tip[0],tip[1],2.6,0,7); ctx.fillStyle="#ff7a5e"; ctx.fill();
      }
    }
  }

  // --- 节点 ---
  if(S.snode||S.slabel){
    for(let i=0;i<nodes.length;++i){
      const p=project(...nodePos(i,true));
      if(S.snode){ ctx.beginPath(); ctx.arc(p[0],p[1],1.8,0,7); ctx.fillStyle="#cfe0f5"; ctx.fill(); }
      if(S.slabel){ ctx.fillStyle="#8b96a8"; ctx.font="10px sans-serif";
        ctx.fillText(String(i),p[0]+3,p[1]-3); }
    }
  }

  // --- 坐标轴 ---
  if(S.sax){
    const o=project(BB.min[0],BB.min[1],BB.min[2]);
    const L=BB.size*0.12;
    const ax=[[[L,0,0],"#ff6b6b","X"],[[0,L,0],"#51cf66","Y"],[[0,0,L],"#4aa3ff","Z"]];
    ctx.font="11px sans-serif";
    for(const [d,c,nm] of ax){
      const q=project(BB.min[0]+d[0],BB.min[1]+d[1],BB.min[2]+d[2]);
      ctx.beginPath(); ctx.moveTo(o[0],o[1]); ctx.lineTo(q[0],q[1]);
      ctx.strokeStyle=c; ctx.lineWidth=2; ctx.stroke();
      ctx.fillStyle=c; ctx.fillText(nm,q[0]+4,q[1]-2);
    }
  }
  updateLegend(f,fmin,fmax);
  updateStatus();
}
function lerpNode(i,j,t){
  const a=nodePos(i,true), b=nodePos(j,true);
  return [a[0]+(b[0]-a[0])*t, a[1]+(b[1]-a[1])*t, a[2]+(b[2]-a[2])*t];
}

// ---- 图例 ----
function updateLegend(f,mn,mx){
  const lc=document.getElementById("lg"), lx=lc.getContext("2d");
  for(let x=0;x<lc.width;++x){ lx.fillStyle=rgb(jet(x/(lc.width-1))); lx.fillRect(x,0,1,lc.height); }
  document.getElementById("lgmin").textContent=F(mn,3);
  document.getElementById("lgmax").textContent=F(mx,3);
  document.getElementById("lgname").textContent=f?(f.name+(f.unit&&f.unit!=="—"?" ("+f.unit+")":"")):"单色";
  document.getElementById("crange").textContent=f?(F(mn,3)+" ~ "+F(mx,3)):"—";
}
function updateStatus(){
  document.getElementById("status").textContent =
    `节点 ${nodes.length} · 杆件 ${beams.length} · 壳 ${shells.length} · 放大 ${S.dispScale.toFixed(0)}×`;
}

// ---- 选中构件 ----
function pick(mx,my){
  let best=-1,bd=8;
  for(let bi=0;bi<beams.length;++bi){
    const b=beams[bi];
    const a=project(...nodePos(b[BNI],true)), c=project(...nodePos(b[BNJ],true));
    const d=segDist(mx,my,a[0],a[1],c[0],c[1]);
    if(d<bd){bd=d;best=bi;}
  }
  S.sel=best;
  const rp=document.getElementById("right");
  if(best>=0){ rp.classList.add("on"); showMember(best); } else rp.classList.remove("on");
  draw();
}
function segDist(px,py,x1,y1,x2,y2){
  const dx=x2-x1,dy=y2-y1; const L2=dx*dx+dy*dy;
  let t=L2>0?((px-x1)*dx+(py-y1)*dy)/L2:0; t=clamp(t,0,1);
  return Math.hypot(px-(x1+t*dx),py-(y1+t*dy));
}
function showMember(i){
  const b=beams[i];
  const rows=[["编号","#"+b[BID]],["类别",KINDNAME[b[BKIND]]||"其它"],
    ["长度 (m)",F(b[BL],3)],["自重 (kN)",F(b[BW],3)],
    ["最大正应力 (kPa)",F(b[BSMAX],2)],["最小正应力 (kPa)",F(b[BSMIN],2)],
    ["von Mises (kPa)",F(b[BVM],2)],["应力比",F(b[BUTIL],3)],
    ["最大弯矩 (kN·m)",F(b[BMMAX],3)],["最大剪力 (kN)",F(b[BVMAX],3)],
    ["轴力 最大/最小 (kN)",F(b[BNMAX],3)+" / "+F(b[BNMIN],3)]];
  document.getElementById("minfo").innerHTML=
    rows.map(r=>`<div class="kv"><span>${r[0]}</span><span>${r[1]}</span></div>`).join("");
  drawDiagram(b);
  // 沿杆极值（重新按 station 扫一遍，避免与 C++ 采样冲突）
  const st=b[BST]||[];
  let Mmax=0,Vmax=0,Nmax=-1e30,Nmin=1e30;
  for(const s of st){ Mmax=Math.max(Mmax,Math.abs(s[SMY]),Math.abs(s[SMZ]));
    Vmax=Math.max(Vmax,Math.abs(s[SVY]),Math.abs(s[SVZ]));
    Nmax=Math.max(Nmax,s[SN]); Nmin=Math.min(Nmin,s[SN]); }
  document.getElementById("mstat").innerHTML=
    [["|M|max (kN·m)",F(Mmax,3)],["|V|max (kN)",F(Vmax,3)],
     ["Nmax (kN)",F(Nmax,3)],["Nmin (kN)",F(Nmin,3)],["分段数",String(st.length)]]
    .map(r=>`<div class="kv"><span>${r[0]}</span><span>${r[1]}</span></div>`).join("");
}
function drawDiagram(b){
  const c=document.getElementById("dg"), g=c.getContext("2d");
  g.clearRect(0,0,c.width,c.height);
  const st=b[BST]||[];
  if(st.length<2){ g.fillStyle="#8b96a8"; g.font="12px sans-serif";
    g.fillText("无分段数据（导出时已关闭）",10,20); return; }
  const H=c.height/3-6, pad=26, W=c.width-pad-8;
  const series=[["N (kN)",s=>s[SN],"#4aa3ff"],
                ["V (kN)",s=>s[SVZ],"#37c76a"],
                ["M (kN·m)",s=>s[SMZ],"#ffb648"]];
  series.forEach((ser,si)=>{
    const top=si*(H+6)+4;
    let mn=1e30,mx=-1e30;
    for(const s of st){ const v=ser[1](s); if(v<mn)mn=v; if(v>mx)mx=v; }
    const sp=(mx-mn)||1;
    // 基线
    const y0=top+H-(0-mn)/sp*H;
    g.strokeStyle="#3a4554"; g.lineWidth=1;
    g.beginPath(); g.moveTo(pad,y0); g.lineTo(pad+W,y0); g.stroke();
    // 曲线（弯矩画在受拉侧：工程习惯把 M 画在受拉边，这里直接按值画）
    g.beginPath();
    st.forEach((s,k)=>{ const x=pad+W*s[SXI], y=top+H-(ser[1](s)-mn)/sp*H;
      k?g.lineTo(x,y):g.moveTo(x,y); });
    g.strokeStyle=ser[2]; g.lineWidth=1.8; g.stroke();
    // 填充
    g.lineTo(pad+W,y0); g.lineTo(pad,y0); g.closePath();
    g.fillStyle=ser[2]+"33"; g.fill();
    g.fillStyle="#dfe4ec"; g.font="11px sans-serif";
    g.fillText(ser[0],4,top+11);
    g.fillStyle="#8b96a8"; g.font="10px sans-serif";
    g.fillText(F(mx,2),4,top+H-1); g.fillText(F(mn,2),4,top+22);
  });
}

// ---- 交互 ----
let dragging=0,lastX=0,lastY=0;
cv.addEventListener("mousedown",e=>{dragging=(e.button===2)?2:1;lastX=e.clientX;lastY=e.clientY;
  cv.classList.add("drag");});
window.addEventListener("mouseup",()=>{dragging=0;cv.classList.remove("drag");});
window.addEventListener("mousemove",e=>{
  if(!dragging)return;
  const dx=e.clientX-lastX, dy=e.clientY-lastY;
  lastX=e.clientX; lastY=e.clientY;
  if(dragging===1){ S.yaw+=dx*0.008; S.pitch=clamp(S.pitch+dy*0.008,-1.5,1.5); }
  else { S.panX+=dx; S.panY+=dy; }
  draw();
});
cv.addEventListener("contextmenu",e=>e.preventDefault());
cv.addEventListener("wheel",e=>{e.preventDefault();
  S.zoom*=Math.exp(-e.deltaY*0.0012); S.zoom=clamp(S.zoom,0.05,60);
  VP.s=(Math.min(cv.getBoundingClientRect().width,cv.getBoundingClientRect().height)/BB.size)*0.72*S.zoom;
  draw();},{passive:false});
cv.addEventListener("click",e=>{
  const r=cv.getBoundingClientRect();
  pick(e.clientX-r.left,e.clientY-r.top);
});
cv.addEventListener("mousemove",e=>{
  if(dragging)return;
  const r=cv.getBoundingClientRect();
  const mx=e.clientX-r.left,my=e.clientY-r.top;
  let txt="";
  for(let i=0;i<nodes.length;++i){
    const p=project(...nodePos(i,true));
    if(Math.hypot(p[0]-mx,p[1]-my)<5){ txt=`节点 ${i}\n位移 (${F(nodes[i][NUX],5)}, ${F(nodes[i][NUY],5)}, ${F(nodes[i][NUZ],5)}) m`; break; }
  }
  const el=document.getElementById("pick");
  if(txt){ el.style.display="block"; el.textContent=txt;
    el.style.left=(mx+12)+"px"; el.style.top=(my+12)+"px"; }
  else el.style.display="none";
});
const VIEWS={iso:[-0.62,0.42],top:[0,1.5707],front:[0,0],side:[1.5707,0]};
document.querySelectorAll("[data-view]").forEach(b=>b.onclick=()=>{
  const v=VIEWS[b.dataset.view]; S.yaw=v[0]; S.pitch=v[1]; S.zoom=1;
  document.getElementById("fit").onclick(); draw();});
document.getElementById("fit").onclick=()=>{
  S.panX=0;S.panY=0;S.zoom=1;
  const r=cv.getBoundingClientRect();
  VP.s=(Math.min(r.width,r.height)/BB.size)*0.72; draw();};

// ---- 面板 ----
(function initModes(){
  const sel=document.getElementById("cmode");
  let opts=`<option value="">单色（不着色）</option>`;
  opts+=`<option value="vm">梁 von Mises 应力</option>`;
  opts+=`<option value="util">梁 应力比</option>`;
  opts+=`<option value="axial">梁 轴力极值</option>`;
  fields.forEach((f,i)=>{ opts+=`<option value="f:${i}">节点场：${f.name}${f.unit&&f.unit!=="—"?" ("+f.unit+")":""}</option>`; });
  sel.innerHTML=opts;
  // 默认选第一个"位移"相关场，最直观
  const pref=fields.findIndex(f=>/位移|U/.test(f.name));
  sel.value = pref>=0?("f:"+pref):(fields.length?"f:0":"");
  S.cmode=sel.value;
  sel.onchange=()=>{S.cmode=sel.value;draw();};
})();
function bindCheck(id,key){ const el=document.getElementById(id);
  el.checked=S[key]; el.onchange=()=>{S[key]=el.checked;draw();}; }
bindCheck("snode","snode"); bindCheck("sbeam","sbeam"); bindCheck("sshell","sshell");
bindCheck("sreact","sreact"); bindCheck("saxes","sax"); bindCheck("slabel","slabel");
bindCheck("deform","deform"); bindCheck("undef","undef");

const ds=document.getElementById("dscale");
function applyScale(){
  // 对数刻度：0..1000 → 1..10000 倍
  S.dispScale=Math.round(Math.pow(10,ds.value/1000*4));
  document.getElementById("dscaletxt").textContent=S.dispScale+"×";
}
ds.oninput=()=>{applyScale();draw();};
(function autoScale(){
  const md=env.maxDisp||0;
  if(md>1e-12){ const want=BB.size*0.06/md; ds.value=clamp(Math.log10(Math.max(want,1))/4*1000,0,1000); }
  applyScale();
})();
document.getElementById("play").onclick=function(){
  S.anim=!S.anim; this.textContent=S.anim?"⏸ 暂停":"▶ 动画";
  this.classList.toggle("act",S.anim);
  if(S.anim){ let t0=performance.now();
    const step=(t)=>{ if(!S.anim)return;
      S.animT=(Math.sin((t-t0)/900)+1)/2*1.0;   // 0..1 往复
      draw(); requestAnimationFrame(step); };
    requestAnimationFrame(step);
  } else { S.animT=1; draw(); }
};
document.getElementById("reset").onclick=()=>{S.anim=false;S.animT=1;
  document.getElementById("play").textContent="▶ 动画";
  document.getElementById("play").classList.remove("act");
  S.yaw=-0.62;S.pitch=0.42;S.zoom=1;S.panX=0;S.panY=0;
  document.getElementById("fit").onclick();};
document.getElementById("close").onclick=()=>{S.sel=-1;
  document.getElementById("right").classList.remove("on");draw();};

// ---- 顶部 / 统计 / 楼层 ----
(function header(){
  document.getElementById("hsub").textContent=DATA.summary||"";
  const ok=document.getElementById("hok");
  ok.textContent=meta.ok?"求解成功":"求解失败"; ok.className="badge "+(meta.ok?"ok":"bad");
  const rs=document.getElementById("hres");
  const r=meta.residual||0;
  rs.textContent="残差 "+r.toExponential(2);
  rs.className="badge "+(r<1e-8?"ok":(r<1e-4?"":"bad"));
  document.getElementById("hsize").textContent=
    `自由度 ${meta.ndof} · nnz ${meta.nnz} · ${(meta.seconds||0).toFixed(3)} s`;
})();
(function stats(){
  const kv=[["节点数",nodes.length],["杆件数",beams.length],["壳单元数",shells.length],
    ["自由度数",meta.ndof],["约束自由度",meta.nconstrained],
    ["总重量 (kN)",F(env.totalWeight,2)],["基底剪力 (kN)",F(env.baseShear,3)],
    ["剪重比",F(env.shearWeightRatio,5)],
    ["最大位移 (m)",F(env.maxDisp,6)],["最大水平位移 (m)",F(env.maxUxy,6)],
    ["最大竖向位移 (m)",F(env.maxUz,6)],
    ["最大层间位移角","1/"+(env.maxDrift?Math.round(1/env.maxDrift):"—")],
    ["最大 von Mises (kPa)",F(env.maxVonMises,2)],
    ["最大支座反力 (kN)",F(env.maxReaction,3)],
    ["最大竖向反力 (kN)",F(env.maxReactionZ,3)]];
  document.getElementById("stats").innerHTML=
    kv.map(r=>`<div class="kv"><span>${r[0]}</span><span>${r[1]}</span></div>`).join("");
})();
(function storyTable(){
  if(!stories.length){ document.getElementById("stories").innerHTML=
    '<div style="color:var(--dim);font-size:11px">无楼层信息</div>'; return; }
  let h='<table><tr><th>层</th><th>标高</th><th>层高</th><th>最大Δ</th><th>位移角</th>'
      +'<th>位移比</th><th>剪力</th><th>剪重比</th></tr>';
  for(const s of stories){
    h+=`<tr><td>${s.story}</td><td>${F(s.z,2)}</td><td>${F(s.h,2)}</td>`
      +`<td>${(s.maxUxy*1000).toFixed(3)}</td><td>1/${s.drift?Math.round(1/s.drift):"—"}</td>`
      +`<td>${F(s.driftRatio,3)}</td><td>${F(s.shear,2)}</td><td>${F(s.shearWeightRatio,5)}</td></tr>`;
  }
  document.getElementById("stories").innerHTML=h+"</table>";
})();
(function warns(){
  const w=[];
  if(meta.residual>1e-6) w.push("解残差偏大："+meta.residual.toExponential(2)+" —— 检查约束是否充分");
  if(!meta.ok) w.push("求解未成功："+(meta.message||""));
  if(env.maxDrift>1/550) w.push("最大层间位移角 1/"+Math.round(1/env.maxDrift)+" 超过 1/550（框架限值）");
  const z=document.getElementById("warnzone");
  z.innerHTML = w.length? w.map(t=>`<div class="warnbox">${t}</div>`).join("") : "";
})();

window.addEventListener("resize",resize);
resize();
</script>
</body>
</html>
)YJKRAW";

}  // namespace post
}  // namespace yjk

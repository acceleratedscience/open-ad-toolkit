const __vite__mapDeps=(i,m=__vite__mapDeps,d=(m.f||(m.f=["assets/DataViewer-zQ2r5EOP.js","assets/index-DwYRJoEC.js","assets/index-BpBKVyNv.css","assets/TextBox-BFX_-vlx.js","assets/TextBox-CUZX0UxP.css","assets/DataViewer--1Nvov8V.css","assets/FileBrowser-BFr-95JS.js","assets/FileBrowser-o1RQTRSH.css","assets/JsonViewer-C2ThBJnO.js","assets/JsonViewer.vue_vue_type_script_setup_true_lang-CZUNuMu-.js","assets/BreadCrumbs-CUWBa4am.js","assets/BreadCrumbs-DvCDsjit.css","assets/MolViewer-ealKfnRt.js","assets/MolViewerStore-VkF0_V-U.js","assets/BaseFetching-CtAuO-DI.js","assets/BaseFetching-CiTUXLMo.css","assets/MolRender3D-C-QWcCb6.js","assets/MolRender3D-BOIri_dQ.css","assets/MolViewer-E2DWdpAy.css","assets/MolsetViewer-Da85EEZq.js","assets/MolsetViewer.vue_vue_type_script_setup_true_lang-BNP5W5IC.js","assets/initRDKit-WaX80hiK.js","assets/16-Oa6gYw0C.js","assets/MolsetViewer-DQo9B4Y9.css","assets/PdfViewer-6NVFSMtC.js","assets/TextViewer--bkp-B6RVTPZp.js","assets/TextViewer--bkp-Dlkg3jVw.css","assets/TextViewer-bhxxDgWZ.js","assets/UnknownViewer-AV8PCFyp.js","assets/_templateModule-rd0AuYhp.js"])))=>i.map(i=>d[i]);
import{d as P,u as x,a as R,b as C,r as y,s as g,c as V,w as h,o as L,e as b,f as k,g as F,h as a,i as E,t as I,j as i,k as c,l as M,m as O,n as B,p as S,q as N,F as j,v as q,x as G,_ as $,y as t}from"./index-DwYRJoEC.js";import{u as U}from"./MolViewerStore-VkF0_V-U.js";import{u as z}from"./BaseFetching-CtAuO-DI.js";import{B as H}from"./BreadCrumbs-CUWBa4am.js";import{B as J}from"./BaseFetchingFile-DZ9OR-uR.js";const Q={key:0,class:"error-msg"},K=q("p",null,"File not found",-1),re=P({__name:"ViewerDispatch",setup(W){const T=x(),_=R(),e=C(),f=U(),l=z();let n=null;const p=y(!1),d=y(!1),u=g(null),D=V(()=>!(e.isDir||!e.data||e.moduleName=="MolViewer"||e.moduleName=="MolGrid")),A=V(()=>"/~/"+e.breadCrumbPathArray.slice(1).slice(0,-1).join("/"));v(),h(()=>_.path,()=>{v()}),h(()=>e.fileType,()=>{m(e.moduleName)}),L(e.clear);async function v(){const s=_.path.replace(/(^\/headless)?\/~(\/)?/,""),r=_.query;b(k.getFile(s,r),{loading:p,onError:o=>{e.setForcedLoading(!1),console.error("Error in getFile()",o)},onSuccess:o=>{if(e.setForcedLoading(!1),e.loadItem(o),e.errCode){n!="error"&&m(null),n="error";return}if(["molset","sdf","smi"].includes(e.fileType||"")){const w=o.data;l.setMolset(w),e.ext=="json"?l.setContext("json"):e.ext=="sdf"?l.setContext("sdf-file"):e.ext=="smi"&&l.setContext("smi-file")}else e.fileType=="smol"?f.setMolData(o.data,"smol"):e.fileType&&["mmol","pdb","cif"].includes(e.fileType)&&f.setMolData(o.data,"mmol");n!="file"&&m(e.moduleName),n="file"}})}function m(s){if(d.value=!1,!s){u.value=null;return}u.value=G(()=>$(Object.assign({"../viewers/DataViewer.vue":()=>t(()=>import("./DataViewer-zQ2r5EOP.js"),__vite__mapDeps([0,1,2,3,4,5])),"../viewers/FileBrowser.vue":()=>t(()=>import("./FileBrowser-BFr-95JS.js"),__vite__mapDeps([6,1,2,7])),"../viewers/JsonViewer.vue":()=>t(()=>import("./JsonViewer-C2ThBJnO.js"),__vite__mapDeps([8,9,1,2,10,11,3,4])),"../viewers/MolViewer.vue":()=>t(()=>import("./MolViewer-ealKfnRt.js").then(r=>r.c),__vite__mapDeps([12,1,2,13,14,15,10,11,16,17,18])),"../viewers/MolsetViewer.vue":()=>t(()=>import("./MolsetViewer-Da85EEZq.js"),__vite__mapDeps([19,20,1,2,14,15,13,10,11,12,16,17,18,21,22,23])),"../viewers/PdfViewer.vue":()=>t(()=>import("./PdfViewer-6NVFSMtC.js"),__vite__mapDeps([24,1,2,10,11])),"../viewers/TextViewer--bkp.vue":()=>t(()=>import("./TextViewer--bkp-B6RVTPZp.js"),__vite__mapDeps([25,1,2,10,11,26])),"../viewers/TextViewer.vue":()=>t(()=>import("./TextViewer-bhxxDgWZ.js"),__vite__mapDeps([27,1,2,10,11,3,4])),"../viewers/UnknownViewer.vue":()=>t(()=>import("./UnknownViewer-AV8PCFyp.js"),__vite__mapDeps([28,1,2,10,11])),"../viewers/_templateModule.vue":()=>t(()=>import("./_templateModule-rd0AuYhp.js"),__vite__mapDeps([29,1,2]))}),`../viewers/${s}.vue`,3).catch(r=>{console.log(111,`../viewers/${s}.vue`,r),d.value=!0}))}return(s,r)=>{const o=F("cv-button");return d.value?(a(),E("div",Q,"The requested module '"+I(i(e).moduleName)+"' was not found.",1)):i(e).forcedLoading||p.value&&!i(e).isDir?(a(),c(J,{key:1,text:"Opening file",failText:"Failed to open file"})):u.value?(a(),c(M(u.value),{key:2})):(a(),E(j,{key:3},[D.value?O("",!0):(a(),c(H,{key:0,pathArray:i(e).breadCrumbPathArray},null,8,["pathArray"])),K,B(o,{size:"small",onClick:r[0]||(r[0]=w=>i(T).push(A.value))},{default:S(()=>[N("Exit")]),_:1})],64))}}});export{re as default};

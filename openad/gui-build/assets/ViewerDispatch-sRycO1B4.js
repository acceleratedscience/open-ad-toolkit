const __vite__mapDeps=(i,m=__vite__mapDeps,d=(m.f||(m.f=["assets/DataViewer-C3OtVIkb.js","assets/index-D-rujQAW.js","assets/index-DfkNSHLR.css","assets/TextBox-DYqn3lGo.js","assets/TextBox-CUZX0UxP.css","assets/DataViewer--1Nvov8V.css","assets/FileBrowser-CekbWbUV.js","assets/FileBrowser-CRInM39Z.css","assets/JsonViewer-x_70-D9y.js","assets/JsonViewer.vue_vue_type_script_setup_true_lang-DBJXQg-x.js","assets/BreadCrumbs-6ctqne71.js","assets/BreadCrumbs-DvCDsjit.css","assets/MolViewer-DMExLisu.js","assets/MolGridStore-ikTBQjQU.js","assets/BaseFetching-DWxuFDXw.js","assets/BaseFetching-DS2kHIcg.css","assets/MolViewer-0Y1SC3um.css","assets/MolsetViewer-DKX261Hp.js","assets/MolsetViewer.vue_vue_type_script_setup_true_lang-DZuB-rAv.js","assets/initRDKit-WaX80hiK.js","assets/16-CAW1pidv.js","assets/MolsetViewer-CWCMe-Q_.css","assets/PdfViewer-CFWLekLx.js","assets/TextViewer--bkp-BS2-ufEn.js","assets/TextViewer--bkp-Dlkg3jVw.css","assets/TextViewer-C6w26QQm.js","assets/UnknownViewer-DqtthTBQ.js","assets/_templateModule-DMoqtx1N.js"])))=>i.map(i=>d[i]);
import{d as P,u as x,a as R,b as C,c as k,r as w,s as b,e as V,w as h,o as I,f as g,g as B,h as F,i as a,j as y,t as L,k as _,l as f,m as M,n as O,p as S,q as N,v as j,F as q,x as G,y as z,_ as U,z as t}from"./index-D-rujQAW.js";import{u as $}from"./MolGridStore-ikTBQjQU.js";import{B as H}from"./BreadCrumbs-6ctqne71.js";import{B as J}from"./BaseFetching-DWxuFDXw.js";const Q={key:0,class:"error-msg"},K=G("p",null,"File not found",-1),te=P({__name:"ViewerDispatch",setup(W){const E=x(),m=R(),e=C(),T=k(),i=$();let l=null;const p=w(!1),c=w(!1),n=b(null),A=V(()=>!(e.isDir||!e.data||e.moduleName=="MolViewer"||e.moduleName=="MolGrid")),D=V(()=>"/~/"+e.breadCrumbPathArray.slice(1).slice(0,-1).join("/"));v(),h(()=>m.path,()=>{v()}),h(()=>e.fileType,(o,r)=>{d(e.moduleName)}),I(e.clear);async function v(){const o=m.path.replace(/(^\/headless)?\/~(\/)?/,""),r=m.query;g(B.getFile(o,r),{loading:p,onError:s=>{console.log("Error in getFile()",s)},onSuccess:s=>{if(e.loadItem(s),e.errCode){l!="error"&&d(null),l="error";return}if(["molset","sdf","smi"].includes(e.fileType||"")){const u=s.data;i.setMolset(u),e.ext=="json"?i.setContext("json"):e.ext=="sdf"?i.setContext("sdf-file"):e.ext=="smi"&&i.setContext("smi-file")}else if(e.fileType=="mol"){const u=s.data;T.setMolData(u)}l!="file"&&d(e.moduleName),l="file"}})}function d(o){if(c.value=!1,!o){n.value=null;return}n.value=z(()=>U(Object.assign({"../viewers/DataViewer.vue":()=>t(()=>import("./DataViewer-C3OtVIkb.js"),__vite__mapDeps([0,1,2,3,4,5])),"../viewers/FileBrowser.vue":()=>t(()=>import("./FileBrowser-CekbWbUV.js"),__vite__mapDeps([6,1,2,7])),"../viewers/JsonViewer.vue":()=>t(()=>import("./JsonViewer-x_70-D9y.js"),__vite__mapDeps([8,9,1,2,10,11,3,4])),"../viewers/MolViewer.vue":()=>t(()=>import("./MolViewer-DMExLisu.js").then(r=>r.c),__vite__mapDeps([12,1,2,13,10,11,14,15,16])),"../viewers/MolsetViewer.vue":()=>t(()=>import("./MolsetViewer-DKX261Hp.js"),__vite__mapDeps([17,18,1,2,13,10,11,12,14,15,16,19,20,21])),"../viewers/PdfViewer.vue":()=>t(()=>import("./PdfViewer-CFWLekLx.js"),__vite__mapDeps([22,1,2,10,11])),"../viewers/TextViewer--bkp.vue":()=>t(()=>import("./TextViewer--bkp-BS2-ufEn.js"),__vite__mapDeps([23,1,2,10,11,24])),"../viewers/TextViewer.vue":()=>t(()=>import("./TextViewer-C6w26QQm.js"),__vite__mapDeps([25,1,2,10,11,3,4])),"../viewers/UnknownViewer.vue":()=>t(()=>import("./UnknownViewer-DqtthTBQ.js"),__vite__mapDeps([26,1,2,10,11])),"../viewers/_templateModule.vue":()=>t(()=>import("./_templateModule-DMoqtx1N.js"),__vite__mapDeps([27,1,2]))}),`../viewers/${o}.vue`,3).catch(()=>{c.value=!0}))}return(o,r)=>{const s=F("cv-button");return c.value?(a(),y("div",Q,"The requested module '"+L(_(e).moduleName)+"' was not found.",1)):p.value&&!_(e).isDir?(a(),f(J,{key:1,text:"Fetching file",failText:"Failed to fetch file"})):n.value?(a(),f(M(n.value),{key:2})):(a(),y(q,{key:3},[A.value?O("",!0):(a(),f(H,{key:0,pathArray:_(e).breadCrumbPathArray},null,8,["pathArray"])),K,S(s,{size:"small",onClick:r[0]||(r[0]=u=>_(E).push(D.value))},{default:N(()=>[j("Exit")]),_:1})],64))}}});export{te as default};

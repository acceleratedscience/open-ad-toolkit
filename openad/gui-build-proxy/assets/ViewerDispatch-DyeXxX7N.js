const __vite__mapDeps=(i,m=__vite__mapDeps,d=(m.f||(m.f=["assets/DataViewer-IaIRZV9F.js","assets/index-M2pymOVI.js","assets/index-Bv2_Ke9_.css","assets/TextBox-Dt-kVu0X.js","assets/TextBox-CUZX0UxP.css","assets/DataViewer--1Nvov8V.css","assets/FileBrowser-Bf52lquy.js","assets/FileBrowser-CRInM39Z.css","assets/JsonViewer-BJ4GJc49.js","assets/JsonViewer.vue_vue_type_script_setup_true_lang-5lDjPVoo.js","assets/BreadCrumbs-BJsAHiin.js","assets/BreadCrumbs-DvCDsjit.css","assets/MolViewer-CKUeXNVw.js","assets/MolGridStore-B_qX5m1O.js","assets/BaseFetching-a6W39JZz.js","assets/BaseFetching-DS2kHIcg.css","assets/MolViewer-0Y1SC3um.css","assets/MolsetViewer-ChXTVrtO.js","assets/MolsetViewer.vue_vue_type_script_setup_true_lang-CGPngmF8.js","assets/initRDKit-WaX80hiK.js","assets/16-B-i9dxBf.js","assets/MolsetViewer-CWCMe-Q_.css","assets/PdfViewer-KmxR1FmO.js","assets/TextViewer--bkp-DPxYcqRs.js","assets/TextViewer--bkp-Dlkg3jVw.css","assets/TextViewer-BplwlV8y.js","assets/UnknownViewer-Cr9zkCEb.js","assets/_templateModule-BWfcyBvX.js"])))=>i.map(i=>d[i]);
import{d as P,u as x,a as R,b as C,c as k,r as w,s as b,e as V,w as h,o as I,f as g,g as B,h as F,i as a,j as y,t as L,k as _,l as f,m as M,n as O,p as S,q as N,v as j,F as q,x as G,y as z,_ as U,z as t}from"./index-M2pymOVI.js";import{u as $}from"./MolGridStore-B_qX5m1O.js";import{B as H}from"./BreadCrumbs-BJsAHiin.js";import{B as J}from"./BaseFetching-a6W39JZz.js";const Q={key:0,class:"error-msg"},K=G("p",null,"File not found",-1),te=P({__name:"ViewerDispatch",setup(W){const E=x(),m=R(),e=C(),T=k(),i=$();let l=null;const p=w(!1),c=w(!1),n=b(null),A=V(()=>!(e.isDir||!e.data||e.moduleName=="MolViewer"||e.moduleName=="MolGrid")),D=V(()=>"/~/"+e.breadCrumbPathArray.slice(1).slice(0,-1).join("/"));v(),h(()=>m.path,()=>{v()}),h(()=>e.fileType,(o,r)=>{d(e.moduleName)}),I(e.clear);async function v(){const o=m.path.replace(/(^\/headless)?\/~(\/)?/,""),r=m.query;g(B.getFile(o,r),{loading:p,onError:s=>{console.log("Error in getFile()",s)},onSuccess:s=>{if(e.loadItem(s),e.errCode){l!="error"&&d(null),l="error";return}if(["molset","sdf","smi"].includes(e.fileType||"")){const u=s.data;i.setMolset(u),e.ext=="json"?i.setContext("json"):e.ext=="sdf"?i.setContext("sdf-file"):e.ext=="smi"&&i.setContext("smi-file")}else if(e.fileType=="mol"){const u=s.data;T.setMolData(u)}l!="file"&&d(e.moduleName),l="file"}})}function d(o){if(c.value=!1,!o){n.value=null;return}n.value=z(()=>U(Object.assign({"../viewers/DataViewer.vue":()=>t(()=>import("./DataViewer-IaIRZV9F.js"),__vite__mapDeps([0,1,2,3,4,5])),"../viewers/FileBrowser.vue":()=>t(()=>import("./FileBrowser-Bf52lquy.js"),__vite__mapDeps([6,1,2,7])),"../viewers/JsonViewer.vue":()=>t(()=>import("./JsonViewer-BJ4GJc49.js"),__vite__mapDeps([8,9,1,2,10,11,3,4])),"../viewers/MolViewer.vue":()=>t(()=>import("./MolViewer-CKUeXNVw.js").then(r=>r.c),__vite__mapDeps([12,1,2,13,10,11,14,15,16])),"../viewers/MolsetViewer.vue":()=>t(()=>import("./MolsetViewer-ChXTVrtO.js"),__vite__mapDeps([17,18,1,2,13,10,11,12,14,15,16,19,20,21])),"../viewers/PdfViewer.vue":()=>t(()=>import("./PdfViewer-KmxR1FmO.js"),__vite__mapDeps([22,1,2,10,11])),"../viewers/TextViewer--bkp.vue":()=>t(()=>import("./TextViewer--bkp-DPxYcqRs.js"),__vite__mapDeps([23,1,2,10,11,24])),"../viewers/TextViewer.vue":()=>t(()=>import("./TextViewer-BplwlV8y.js"),__vite__mapDeps([25,1,2,10,11,3,4])),"../viewers/UnknownViewer.vue":()=>t(()=>import("./UnknownViewer-Cr9zkCEb.js"),__vite__mapDeps([26,1,2,10,11])),"../viewers/_templateModule.vue":()=>t(()=>import("./_templateModule-BWfcyBvX.js"),__vite__mapDeps([27,1,2]))}),`../viewers/${o}.vue`,3).catch(()=>{c.value=!0}))}return(o,r)=>{const s=F("cv-button");return c.value?(a(),y("div",Q,"The requested module '"+L(_(e).moduleName)+"' was not found.",1)):p.value&&!_(e).isDir?(a(),f(J,{key:1,text:"Fetching file",failText:"Failed to fetch file"})):n.value?(a(),f(M(n.value),{key:2})):(a(),y(q,{key:3},[A.value?O("",!0):(a(),f(H,{key:0,pathArray:_(e).breadCrumbPathArray},null,8,["pathArray"])),K,S(s,{size:"small",onClick:r[0]||(r[0]=u=>_(E).push(D.value))},{default:N(()=>[j("Exit")]),_:1})],64))}}});export{te as default};
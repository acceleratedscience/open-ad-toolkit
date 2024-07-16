import{d as P,u as x,a as R,b as C,c as k,r as w,s as b,e as V,w as h,o as I,f as g,g as B,h as F,i as a,j as y,t as L,k as _,l as f,m as M,n as O,p as S,q as N,v as j,F as q,x as G,y as z,_ as U,z as t}from"./index-juGIK60X.js";import{u as $}from"./MolGridStore-YRzXPe7r.js";import{B as H}from"./BreadCrumbs-5PNbBUT3.js";import{B as J}from"./BaseFetching-muE6g1O5.js";const Q={key:0,class:"error-msg"},K=G("p",null,"File not found",-1),te=P({__name:"ViewerDispatch",setup(W){const E=x(),m=R(),e=C(),T=k(),i=$();let l=null;const p=w(!1),c=w(!1),n=b(null),A=V(()=>!(e.isDir||!e.data||e.moduleName=="MolViewer"||e.moduleName=="MolGrid")),D=V(()=>"/~/"+e.breadCrumbPathArray.slice(1).slice(0,-1).join("/"));v(),h(()=>m.path,()=>{v()}),h(()=>e.fileType,(o,r)=>{d(e.moduleName)}),I(e.clear);async function v(){const o=m.path.replace(/(^\/headless)?\/~(\/)?/,""),r=m.query;g(B.getFile(o,r),{loading:p,onError:s=>{console.log("Error in getFile()",s)},onSuccess:s=>{if(e.loadItem(s),e.errCode){l!="error"&&d(null),l="error";return}if(["molset","sdf","smi"].includes(e.fileType||"")){const u=s.data;i.setMolset(u),e.ext=="json"?i.setContext("json"):e.ext=="sdf"?i.setContext("sdf-file"):e.ext=="smi"&&i.setContext("smi-file")}else if(e.fileType=="mol"){const u=s.data;T.setMolData(u)}l!="file"&&d(e.moduleName),l="file"}})}function d(o){if(c.value=!1,!o){n.value=null;return}n.value=z(()=>U(Object.assign({"../viewers/DataViewer.vue":()=>t(()=>import("./DataViewer-cSqhEnNK.js"),__vite__mapDeps([0,1,2,3,4,5])),"../viewers/FileBrowser.vue":()=>t(()=>import("./FileBrowser--NSmbFk2.js"),__vite__mapDeps([6,1,2,7])),"../viewers/JsonViewer.vue":()=>t(()=>import("./JsonViewer-rVWZIkjD.js"),__vite__mapDeps([8,9,1,2,10,11,3,4])),"../viewers/MolViewer.vue":()=>t(()=>import("./MolViewer-8oo8dSmH.js").then(r=>r.c),__vite__mapDeps([12,1,2,13,10,11,14,15,16])),"../viewers/MolsetViewer.vue":()=>t(()=>import("./MolsetViewer-sGE6aZrH.js"),__vite__mapDeps([17,18,1,2,13,10,11,12,14,15,16,19,20,21])),"../viewers/PdfViewer.vue":()=>t(()=>import("./PdfViewer-q7bjtfPQ.js"),__vite__mapDeps([22,1,2,10,11])),"../viewers/TextViewer--bkp.vue":()=>t(()=>import("./TextViewer--bkp-btuqDWqW.js"),__vite__mapDeps([23,1,2,10,11,24])),"../viewers/TextViewer.vue":()=>t(()=>import("./TextViewer-j6hPGBGz.js"),__vite__mapDeps([25,1,2,10,11,3,4])),"../viewers/UnknownViewer.vue":()=>t(()=>import("./UnknownViewer-4UOZT9e7.js"),__vite__mapDeps([26,1,2,10,11])),"../viewers/_templateModule.vue":()=>t(()=>import("./_templateModule-KlGO6hAT.js"),__vite__mapDeps([27,1,2]))}),`../viewers/${o}.vue`).catch(()=>{c.value=!0}))}return(o,r)=>{const s=F("cv-button");return c.value?(a(),y("div",Q,"The requested module '"+L(_(e).moduleName)+"' was not found.",1)):p.value&&!_(e).isDir?(a(),f(J,{key:1,text:"Fetching file",failText:"Failed to fetch file"})):n.value?(a(),f(M(n.value),{key:2})):(a(),y(q,{key:3},[A.value?O("",!0):(a(),f(H,{key:0,pathArray:_(e).breadCrumbPathArray},null,8,["pathArray"])),K,S(s,{size:"small",onClick:r[0]||(r[0]=u=>_(E).push(D.value))},{default:N(()=>[j("Exit")]),_:1})],64))}}});export{te as default};
function __vite__mapDeps(indexes) {
  if (!__vite__mapDeps.viteFileDeps) {
    __vite__mapDeps.viteFileDeps = ["assets/DataViewer-cSqhEnNK.js","assets/index-juGIK60X.js","assets/index-FiPV9bQR.css","assets/TextBox-096FS_mU.js","assets/TextBox-lGV9FMT8.css","assets/DataViewer-PtTb6L_F.css","assets/FileBrowser--NSmbFk2.js","assets/FileBrowser-kSJzN_WV.css","assets/JsonViewer-rVWZIkjD.js","assets/JsonViewer.vue_vue_type_script_setup_true_lang-4fSEEgnR.js","assets/BreadCrumbs-5PNbBUT3.js","assets/BreadCrumbs-7wg7I4rU.css","assets/MolViewer-8oo8dSmH.js","assets/MolGridStore-YRzXPe7r.js","assets/BaseFetching-muE6g1O5.js","assets/BaseFetching-0tpByHIP.css","assets/MolViewer-NGNUgt7p.css","assets/MolsetViewer-sGE6aZrH.js","assets/MolsetViewer.vue_vue_type_script_setup_true_lang-xyjV8fpT.js","assets/initRDKit-FDLbkFFb.js","assets/16-e0gwu24e.js","assets/MolsetViewer-7lUbl9gg.css","assets/PdfViewer-q7bjtfPQ.js","assets/TextViewer--bkp-btuqDWqW.js","assets/TextViewer--bkp-5ZIN41cF.css","assets/TextViewer-j6hPGBGz.js","assets/UnknownViewer-4UOZT9e7.js","assets/_templateModule-KlGO6hAT.js"]
  }
  return indexes.map((i) => __vite__mapDeps.viteFileDeps[i])
}

import{d as g,a as _,r as s,I as h,e as f,A as k,h as t,k as l,i,t as m,m as v,F as y,v as B}from"./index-CPNoIySf.js";import{u as M,B as S}from"./BaseFetching-BkwSakDF.js";import{_ as C}from"./MolsetViewer.vue_vue_type_script_setup_true_lang-qNmCx1_v.js";import"./MolViewerStore-Daj2yBYh.js";import"./BreadCrumbs-DvHG2Qbw.js";import"./MolViewer-9zCwAbA7.js";import"./MolRender3D-B-xg-nr8.js";import"./initRDKit-WaX80hiK.js";import"./16-Gy2UeNuS.js";const E=B("div",{class:"error-msg"},"Something went wrong loading this molecule set.",-1),F={key:0,class:"status-msg"},P=g({__name:"MolsetPage",props:{cacheId:{}},setup(u){const p=_(),a=M(),r=s(!1),e=s(""),n=s(null),d=u;return h(()=>{const c=p.query;f(k.getMolset(+d.cacheId,c),{onSuccess:o=>{a.setMolset(o),a.setContext(null)},onError:o=>{console.log("Error in getMolset()",o)},loading:r,loadingError:e,status:n})}),(c,o)=>r.value?(t(),l(S,{key:0})):e.value?(t(),i(y,{key:1},[E,e.value?(t(),i("div",F,m(n.value)+": "+m(e.value),1)):v("",!0)],64)):(t(),l(C,{key:2,retainCache:!0}))}});export{P as default};

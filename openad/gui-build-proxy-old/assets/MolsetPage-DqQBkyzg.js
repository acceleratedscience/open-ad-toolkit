import{d as g,a as _,r as s,G as f,f as h,A as k,i as t,l,j as i,t as m,n as v,F as y,x as B}from"./index-M2pymOVI.js";import{u as M}from"./MolGridStore-B_qX5m1O.js";import{_ as S}from"./MolsetViewer.vue_vue_type_script_setup_true_lang-CGPngmF8.js";import{B as x}from"./BaseFetching-a6W39JZz.js";import"./BreadCrumbs-BJsAHiin.js";import"./MolViewer-CKUeXNVw.js";import"./initRDKit-WaX80hiK.js";import"./16-B-i9dxBf.js";const C=B("div",{class:"error-msg"},"Something went wrong loading this molecule set.",-1),E={key:0,class:"status-msg"},j=g({__name:"MolsetPage",props:{cacheId:{}},setup(u){const p=_(),a=M(),r=s(!1),e=s(""),n=s(null),d=u;return f(()=>{const c=p.query;h(k.getMolset(+d.cacheId,c),{onSuccess:o=>{a.setMolset(o),a.setContext(null)},onError:o=>{console.log("Error in getMolset()",o)},loading:r,loadingError:e,status:n})}),(c,o)=>r.value?(t(),l(x,{key:0})):e.value?(t(),i(y,{key:1},[C,e.value?(t(),i("div",E,m(n.value)+": "+m(e.value),1)):v("",!0)],64)):(t(),l(S,{key:2,retainCache:!0}))}});export{j as default};

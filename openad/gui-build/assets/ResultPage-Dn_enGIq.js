import{d as g,a as f,r as o,H as k,e as v,J as C,h as t,k as _,i as n,t as a,m as x,F as d,v as s,q as y}from"./index-jyqcHVzo.js";import{u as b,B}from"./BaseFetching-Bsm3mKXj.js";import{_ as T}from"./MolsetViewer.vue_vue_type_script_setup_true_lang-tg1hoiCb.js";import"./MolViewerStore-BLw9AN1l.js";import"./BreadCrumbs-DPqFwJMj.js";import"./MolViewer-BcMJI8ag.js";import"./MolRender3D-DyI99qNw.js";import"./initRDKit-WaX80hiK.js";import"./16-KXY1zdth.js";const S=s("h3",null,"Result",-1),F=s("p",null,[y(" This page displays any data result from the CLI."),s("br"),y(" There is currently no result stored in memory. ")],-1),M=s("p",null,"To give it a try, run in your terminal:",-1),R=s("span",{class:"code",style:{"margin-top":"8px",display:"inline-block"}},"set context gt4sd",-1),q=s("br",null,null,-1),w=s("span",{class:"code",style:{"margin-top":"8px",display:"inline-block"}},"search for similar molecules to 'C1(C(=C)C([O-])C1C)=O'",-1),A=s("br",null,null,-1),E=s("span",{class:"code",style:{"margin-top":"8px",display:"inline-block"}},"result open",-1),N={key:2},V=s("div",{class:"error-msg"},"Something went wrong loading this molecule set.",-1),D={key:0,class:"status-msg"},z=g({__name:"ResultPage",setup(G){const h=f(),i=b(),c=o(!1),u=o(!1),l=o(""),p=o(null),r=o("");return k(()=>{const m=h.query;v(C.getResult(m),{onSuccess:e=>{if(e&&e.type=="empty"){c.value=!0;return}else e&&e.type=="molset"?(i.setMolset(e.data),i.setContext("result-mols")):e&&e.type=="data"?r.value=e.data:console.log("Unknown data type:",e)},loading:u,loadingError:l,status:p})}),(m,e)=>u.value?(t(),_(B,{key:0})):c.value?(t(),n(d,{key:1},[S,F,M,R,q,w,A,E],64)):r.value?(t(),n("pre",N,a(r.value),1)):l.value?(t(),n(d,{key:3},[V,l.value?(t(),n("div",D,a(p.value)+": "+a(l.value),1)):x("",!0)],64)):(t(),_(T,{key:4}))}});export{z as default};

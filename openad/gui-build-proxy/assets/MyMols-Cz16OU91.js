import{d as f,r as l,I as v,e as g,A as M,h as e,k as _,i as o,j as S,F as n,v as s,n as b,J as w,q as m,m as y,D as B,E as x,G as I}from"./index-CoGQGb_u.js";import{u as V,B as E}from"./BaseFetching-KBmQDor1.js";import{u as F}from"./MolViewerStore-GYm8pbqh.js";import{_ as N}from"./MolsetViewer.vue_vue_type_script_setup_true_lang-yaaIEjEL.js";import"./BreadCrumbs-ULNMOtem.js";import"./MolViewer-BFbwHF1f.js";import"./MolRender3D--ExueadL.js";import"./initRDKit-WaX80hiK.js";import"./16-DgTnHlVr.js";const r=a=>(B("data-v-4bbd9a41"),a=a(),x(),a),C={key:0,class:"error-msg"},G=r(()=>s("p",null,"You haven't saved any molecules yet.",-1)),T=r(()=>s("p",null,"To add molecules to your working set, click the bookmark icon on a molecule, or run in your terminal:",-1)),q=r(()=>s("span",{class:"code",style:{"margin-top":"8px",display:"inline-block"}},"add molecule <identifier>",-1)),A=r(()=>s("br",null,null,-1)),j={key:1,id:"about-msg"},D=r(()=>s("br",null,null,-1)),J=f({__name:"MyMols",setup(a){const c=V(),k=F(),d=l(!0),i=l(""),h=l(null),u=l(!1);return v(()=>{const p=c._setUrlQuery();g(M.getMolset_mymols(p),{onSuccess:t=>{t=="empty"?u.value=!0:(c.setMolset(t),c.setContext("my-mols"))},onError:t=>{console.log("Error in getMyMols()",t)},loading:d,status:h,loadingError:i})}),(p,t)=>d.value?(e(),_(E,{key:0})):(e(),o(n,{key:1},[S(k).molFromMolset?y("",!0):(e(),o(n,{key:0},[s("h3",null,[b(w,{icon:"icn-bookmark-full"}),m("My Molecules")]),i.value?(e(),o("p",C,"Something went wrong")):(e(),o(n,{key:1},[u.value?(e(),o(n,{key:0},[G,T,q,A],64)):(e(),o("div",j,[m(" This is your working set of molecules, it is cleared at the end of your session."),D,m(" If you want to preserve this molecule set, you can save it to your workspace under actions. ")]))],64))],64)),!u.value&&!i.value?(e(),_(N,{key:1})):y("",!0)],64))}}),P=I(J,[["__scopeId","data-v-4bbd9a41"]]);export{P as default};
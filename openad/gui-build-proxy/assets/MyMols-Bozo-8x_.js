import{d as h,c as v,r as l,G as g,f as M,A as S,i as e,l as _,j as o,k as b,F as n,x as s,p as w,H as x,v as m,n as y,C as B,D as V,E}from"./index-DOCjoid0.js";import{u as F}from"./MolGridStore-Cio4AaND.js";import{_ as I}from"./MolsetViewer.vue_vue_type_script_setup_true_lang-CJ8gucVC.js";import{B as C}from"./BaseFetching-sic-mAZ7.js";import"./BreadCrumbs-CFfvtUBi.js";import"./MolViewer-xb3arG26.js";import"./initRDKit-WaX80hiK.js";import"./16-DkEyaGLy.js";const r=a=>(B("data-v-4bbd9a41"),a=a(),V(),a),N={key:0,class:"error-msg"},G=r(()=>s("p",null,"You haven't saved any molecules yet.",-1)),T=r(()=>s("p",null,"To add molecules to your working set, click the bookmark icon on a molecule, or run in your terminal:",-1)),A=r(()=>s("span",{class:"code",style:{"margin-top":"8px",display:"inline-block"}},"add molecule <identifier>",-1)),j=r(()=>s("br",null,null,-1)),q={key:1,id:"about-msg"},D=r(()=>s("br",null,null,-1)),H=h({__name:"MyMols",setup(a){const c=F(),k=v(),d=l(!0),i=l(""),f=l(null),u=l(!1);return g(()=>{const p=c._setUrlQuery();M(S.getMolset_mymols(p),{onSuccess:t=>{t=="empty"?u.value=!0:(c.setMolset(t),c.setContext("my-mols"))},onError:t=>{console.log("Error in getMyMols()",t)},loading:d,status:f,loadingError:i})}),(p,t)=>d.value?(e(),_(C,{key:0})):(e(),o(n,{key:1},[b(k).molFromMolset?y("",!0):(e(),o(n,{key:0},[s("h3",null,[w(x,{icon:"icn-bookmark-full"}),m("My Molecules")]),i.value?(e(),o("p",N,"Something went wrong")):(e(),o(n,{key:1},[u.value?(e(),o(n,{key:0},[G,T,A,j],64)):(e(),o("div",q,[m(" This is your working set of molecules, it is cleared at the end of your session."),D,m(" If you want to preserve this molecule set, you can save it to your workspace under actions. ")]))],64))],64)),!u.value&&!i.value?(e(),_(I,{key:1})):y("",!0)],64))}}),O=E(H,[["__scopeId","data-v-4bbd9a41"]]);export{O as default};

import{d as I,L as k,c as B,H as D,g as P,W as V,h as s,k as $,p as _,v as o,t as c,j as v,X as g,i as t,F as f,A as y,S as M,Y as p,q as i,C as z,D as F,E}from"./index-B1am-EuC.js";import{u as L}from"./MolViewerStore-BS84-6gn.js";const n=l=>(z("data-v-56c4f25a"),l=l(),F(),l),N=n(()=>o("h3",null,"Identifiers!",-1)),j={class:"param-wrap"},q=["data-copy"],A=n(()=>o("div",{class:"filler"},null,-1)),H={class:"val"},O=n(()=>o("hr",null,null,-1)),T=n(()=>o("h3",null,"Properties",-1)),U={class:"param-wrap"},W=["data-copy"],X=n(()=>o("div",{class:"filler"},null,-1)),Y={class:"val"},G=I({__name:"ModalMolPreview",emits:["mounted"],setup(l,{emit:S}){const u=k(),w=L(),x=S,r=B(()=>u.data);D(()=>{x("mounted")});async function b(){w.setMolFromMolsetIndex(r.value.index),u.hide()}return(J,K)=>{const C=P("cv-modal"),d=V("click-to-copy");return s(),$(C,{visible:v(u).visible,size:"xs",onPrimaryClick:b,class:"scrollable"},{title:_(()=>[o("h2",null,c(v(g)(r.value.identifiers.name)||"Unnamed Molecule"),1)]),content:_(()=>{var m,h;return[N,o("div",j,[(s(!0),t(f,null,y((m=r.value)==null?void 0:m.identifiers,(e,a)=>(s(),t("div",{key:a,class:M({empty:!e})},[p((s(),t("div",{"data-copy":`${a}: ${e}`,class:"key"},[i(c(a),1)],8,q)),[[d,!!e]]),A,p((s(),t("div",H,[i(c(e||"-"),1)])),[[d,!!e]])],2))),128))]),O,T,o("div",U,[(s(!0),t(f,null,y((h=r.value)==null?void 0:h.properties,(e,a)=>(s(),t("div",{key:a,class:M({empty:!e&&e!==0})},[p((s(),t("div",{"data-copy":`${a}: ${e}`,class:"key"},[i(c(a),1)],8,W)),[[d,!!e||e===0]]),X,p((s(),t("div",Y,[i(c(e||e===0?e:"-"),1)])),[[d,!!e||e===0]])],2))),128))])]}),"secondary-button":_(()=>[i("Cancel")]),"primary-button":_(()=>[i("Open")]),_:1},8,["visible"])}}}),Z=E(G,[["__scopeId","data-v-56c4f25a"]]);export{Z as default};

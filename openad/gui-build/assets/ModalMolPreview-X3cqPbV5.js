import{d as C,z as k,c as B,I as z,g as D,W as P,h as s,k as V,p as _,v as o,t as c,j as v,X as $,i as t,F as f,B as y,S as M,Y as p,q as i,D as g,E as F,G as E}from"./index-D46wLRWb.js";import{u as N}from"./MolViewerStore-COiJdULc.js";const n=l=>(g("data-v-56c4f25a"),l=l(),F(),l),j=n(()=>o("h3",null,"Identifiers!",-1)),q={class:"param-wrap"},G=["data-copy"],L=n(()=>o("div",{class:"filler"},null,-1)),O={class:"val"},T=n(()=>o("hr",null,null,-1)),U=n(()=>o("h3",null,"Properties",-1)),W={class:"param-wrap"},X=["data-copy"],Y=n(()=>o("div",{class:"filler"},null,-1)),A={class:"val"},H=C({__name:"ModalMolPreview",emits:["mounted"],setup(l,{emit:S}){const u=k(),w=N(),x=S,r=B(()=>u.data);z(()=>{x("mounted")});async function b(){w.setMolFromMolsetIndex(r.value.index),u.hide()}return(J,K)=>{const I=D("cv-modal"),d=P("click-to-copy");return s(),V(I,{visible:v(u).visible,size:"xs",onPrimaryClick:b,class:"scrollable"},{title:_(()=>[o("h2",null,c(v($)(r.value.identifiers.name)||"Unnamed Molecule"),1)]),content:_(()=>{var m,h;return[j,o("div",q,[(s(!0),t(f,null,y((m=r.value)==null?void 0:m.identifiers,(e,a)=>(s(),t("div",{key:a,class:M({empty:!e})},[p((s(),t("div",{"data-copy":`${a}: ${e}`,class:"key"},[i(c(a),1)],8,G)),[[d,!!e]]),L,p((s(),t("div",O,[i(c(e||"-"),1)])),[[d,!!e]])],2))),128))]),T,U,o("div",W,[(s(!0),t(f,null,y((h=r.value)==null?void 0:h.properties,(e,a)=>(s(),t("div",{key:a,class:M({empty:!e&&e!==0})},[p((s(),t("div",{"data-copy":`${a}: ${e}`,class:"key"},[i(c(a),1)],8,X)),[[d,!!e||e===0]]),Y,p((s(),t("div",A,[i(c(e||e===0?e:"-"),1)])),[[d,!!e||e===0]])],2))),128))])]}),"secondary-button":_(()=>[i("Cancel")]),"primary-button":_(()=>[i("Open")]),_:1},8,["visible"])}}}),Z=E(H,[["__scopeId","data-v-56c4f25a"]]);export{Z as default};

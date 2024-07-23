import{d as g,M as O,b as R,u as j,e as n,$ as w,r as z,G as J,h as p,i as u,l as L,q as o,v as t,t as i,k as r,U,j as v,x as _,n as T,p as b,F as V,N as q}from"./index-SnNRlG3y.js";const E={key:0,class:"error-msg"},$=g({__name:"ModalViewer",emits:["mounted"],setup(G,{emit:F}){const y=O(),e=R(),d=j(),h=F,x=n(()=>e.fileType?w[e.fileType]:null),N=n(()=>e.defaultFileType?w[e.defaultFileType]:null),f=[["TextViewer","text viewer","text"],["JsonViewer","JSON viewer","json"],["DataViewer","data viewer","data"],["MolViewer","molecule viewer","mol"],["MolsetViewer","molecule set viewer","molset"]],l=z(e.moduleName),S=n(()=>e.fileTypeOverride&&l.value==e.moduleName?"Reset":"Switch"),k=n(()=>l.value==e.moduleName);J(()=>h("mounted"));async function D(){if(l.value==e.defaultModuleName)d.push({path:d.currentRoute.value.path});else{const c=f.filter(([m])=>m==l.value),a=c?c[0][2]:null;a&&d.push("?use="+a)}y.hide()}return(c,a)=>{const m=p("cv-dropdown-item"),M=p("cv-dropdown"),B=p("cv-modal");return u(),L(B,{visible:r(y).visible,size:"xs",onPrimaryClick:D,primaryButtonDisabled:k.value},{title:o(()=>[t(i(r(U)(r(e).fileType??""))+" viewer",1)]),content:o(()=>[r(e).fileTypeOverride?(u(),v("p",E,[t(" This is a "),_("b",null,i(N.value),1),t(" file, but you are currently viewing it in the "),_("b",null,i(x.value),1),t(" viewer. ")])):T("",!0),b(M,{modelValue:l.value,"onUpdate:modelValue":a[0]||(a[0]=s=>l.value=s),label:"Select viewer"},{default:o(()=>[(u(),v(V,null,q(f,([s,C])=>b(m,{key:s,value:s},{default:o(()=>[t(i(C)+" ",1),s==r(e).moduleName?(u(),v(V,{key:0},[t("(default)")],64)):T("",!0)]),_:2},1032,["value"])),64))]),_:1},8,["modelValue"])]),"secondary-button":o(()=>[t("Cancel")]),"primary-button":o(()=>[t(i(S.value),1)]),_:1},8,["visible","primaryButtonDisabled"])}}});export{$ as default};

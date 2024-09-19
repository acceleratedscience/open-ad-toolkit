import{d as g,L,b as O,u as R,c as n,a0 as w,r as j,H as z,g as p,h as u,k as J,p as o,q as t,t as i,j as r,X as q,i as v,v as _,m as T,n as b,F as V,A}from"./index-jyqcHVzo.js";const E={key:0,class:"error-msg"},U=g({__name:"ModalViewer",emits:["mounted"],setup(H,{emit:F}){const y=L(),e=O(),h=F,d=R(),x=n(()=>e.fileType?w[e.fileType]:null),S=n(()=>e.defaultFileType?w[e.defaultFileType]:null),f=[["TextViewer","text viewer","text"],["JsonViewer","JSON viewer","json"],["DataViewer","data viewer","data"],["MolViewer","molecule viewer","mol"],["MolsetViewer","molecule set viewer","molset"]],l=j(e.moduleName),k=n(()=>e.fileTypeOverride&&l.value==e.moduleName?"Reset":"Switch"),N=n(()=>l.value==e.moduleName);z(()=>h("mounted"));async function D(){if(l.value==e.defaultModuleName)d.push({path:d.currentRoute.value.path});else{const c=f.filter(([m])=>m==l.value),a=c?c[0][2]:null;a&&d.push("?use="+a)}y.hide()}return(c,a)=>{const m=p("cv-dropdown-item"),B=p("cv-dropdown"),C=p("cv-modal");return u(),J(C,{visible:r(y).visible,size:"xs",onPrimaryClick:D,primaryButtonDisabled:N.value},{title:o(()=>[t(i(r(q)(r(e).fileType??""))+" viewer",1)]),content:o(()=>[r(e).fileTypeOverride?(u(),v("p",E,[t(" This is a "),_("b",null,i(S.value),1),t(" file, but you are currently viewing it in the "),_("b",null,i(x.value),1),t(" viewer. ")])):T("",!0),b(B,{modelValue:l.value,"onUpdate:modelValue":a[0]||(a[0]=s=>l.value=s),label:"Select viewer"},{default:o(()=>[(u(),v(V,null,A(f,([s,M])=>b(m,{key:s,value:s},{default:o(()=>[t(i(M)+" ",1),s==r(e).moduleName?(u(),v(V,{key:0},[t("(default)")],64)):T("",!0)]),_:2},1032,["value"])),64))]),_:1},8,["modelValue"])]),"secondary-button":o(()=>[t("Cancel")]),"primary-button":o(()=>[t(i(k.value),1)]),_:1},8,["visible","primaryButtonDisabled"])}}});export{U as default};

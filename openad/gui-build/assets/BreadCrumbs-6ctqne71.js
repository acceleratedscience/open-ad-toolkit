import{d as w,L as M,b as x,M as $,r as p,e as F,G as N,h as V,i as e,j as t,x as k,k as f,t as i,n,F as b,N as j,l as v,q as _,v as g,B as C,O as B,P as E,g as I,K as L,E as O}from"./index-D-rujQAW.js";const W={key:0},q={key:3},z={class:"bc-right"},D=w({__name:"BreadCrumbs",props:{pathArray:{},fileType:{}},setup(S){M();const y=x(),A=$(),s=S,l=p(null),a=p(!0),u=p(!1),d=F(()=>s.fileType??y.fileType);N(()=>{l.value&&l.value.offsetWidth<l.value.scrollWidth&&(u.value=!0)});function m(){a.value=!a.value}return(T,c)=>{const h=V("router-link");return e(),t("div",{id:"breadcrumbs-wrap",class:B({truncate:a.value})},[k("div",{id:"breadcrumbs",ref_key:"$breadcrumbs",ref:l,class:B({truncate:a.value,"needs-truncated":u.value})},[d.value?(e(),t("button",{key:0,id:"file-type",onClick:c[0]||(c[0]=r=>f(A).display("ModalViewer"))},i(d.value),1)):n("",!0),(e(!0),t(b,null,j(s.pathArray,(r,o)=>(e(),t(b,{key:o},[o==s.pathArray.length-1?(e(),t("span",W,i(r),1)):o===0?(e(),v(h,{key:1,to:"/~/",class:"dumb"},{default:_(()=>[g(i(r),1)]),_:2},1024)):(e(),v(h,{key:2,to:"/~/"+s.pathArray.slice(1,o+1).join("/"),class:"dumb"},{default:_(()=>[g(i(r),1)]),_:2},1032,["to"])),o<s.pathArray.length-1?(e(),t("span",q,"  ›  ")):n("",!0)],64))),128)),u.value&&!a.value?(e(),t("a",{key:1,href:"#",class:"toggle-hide",onClick:C(m,["prevent"])},"hide")):n("",!0)],2),u.value&&a.value?(e(),t("a",{key:0,href:"#",class:"toggle-show",onClick:C(m,["prevent"])},"show")):n("",!0),k("div",z,[E(T.$slots,"default",{},void 0,!0),d.value?(e(),v(L,{key:0,icon:"icn-doc",iconHover:"icn-doc-full",btnStyle:"soft",mini:"",onClick:c[1]||(c[1]=r=>f(I).openFileOS(f(y).pathAbsolute))})):n("",!0)])],2)}}}),H=O(D,[["__scopeId","data-v-8c1c688c"]]);export{H as B};

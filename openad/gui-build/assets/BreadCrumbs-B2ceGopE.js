import{d as w,R as M,b as $,L as x,r as p,c as F,H as V,g as N,h as e,i as t,v as k,j as f,t as i,m as n,F as b,A as j,k as v,p as _,q as g,B as C,S,T as E,f as H,Q as I,E as L}from"./index-EqhxGKMg.js";const W={key:0},q={key:3},z={class:"bc-right"},D=w({__name:"BreadCrumbs",props:{pathArray:{},fileType:{}},setup(B){M();const y=$(),A=x(),s=B,l=p(null),a=p(!0),u=p(!1),d=F(()=>s.fileType??y.fileType);V(()=>{l.value&&l.value.offsetWidth<l.value.scrollWidth&&(u.value=!0)});function m(){a.value=!a.value}return(T,c)=>{const h=N("router-link");return e(),t("div",{id:"breadcrumbs-wrap",class:S({truncate:a.value})},[k("div",{id:"breadcrumbs",ref_key:"$breadcrumbs",ref:l,class:S({truncate:a.value,"needs-truncated":u.value})},[d.value?(e(),t("button",{key:0,id:"file-type",onClick:c[0]||(c[0]=r=>f(A).display("ModalViewer"))},i(d.value),1)):n("",!0),(e(!0),t(b,null,j(s.pathArray,(r,o)=>(e(),t(b,{key:o},[o==s.pathArray.length-1?(e(),t("span",W,i(r),1)):o===0?(e(),v(h,{key:1,to:"/~/",class:"dumb"},{default:_(()=>[g(i(r),1)]),_:2},1024)):(e(),v(h,{key:2,to:"/~/"+s.pathArray.slice(1,o+1).join("/"),class:"dumb"},{default:_(()=>[g(i(r),1)]),_:2},1032,["to"])),o<s.pathArray.length-1?(e(),t("span",q,"  ›  ")):n("",!0)],64))),128)),u.value&&!a.value?(e(),t("a",{key:1,href:"#",class:"toggle-hide",onClick:C(m,["prevent"])},"hide")):n("",!0)],2),u.value&&a.value?(e(),t("a",{key:0,href:"#",class:"toggle-show",onClick:C(m,["prevent"])},"show")):n("",!0),k("div",z,[E(T.$slots,"default",{},void 0,!0),d.value?(e(),v(I,{key:0,icon:"icn-doc",iconHover:"icn-doc-full",btnStyle:"soft",mini:"",onClick:c[1]||(c[1]=r=>f(H).openFileOS(f(y).pathAbsolute))})):n("",!0)])],2)}}}),Q=L(D,[["__scopeId","data-v-8c1c688c"]]);export{Q as B};
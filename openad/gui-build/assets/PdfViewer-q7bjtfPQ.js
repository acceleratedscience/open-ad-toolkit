import{d as u,u as p,a as c,b as d,g as l,h as f,i as _,j as h,p as n,q as r,I as m,k as t,v as s,F as b,x as o}from"./index-juGIK60X.js";import{B as C}from"./BreadCrumbs-5PNbBUT3.js";const F=o("h3",null,"PDF files are not yet supported.",-1),k=o("p",null,[s("This file should have opened in your default PDF viewer instead."),o("br"),s("If nothing happens, click the button below.")],-1),v=o("br",null,null,-1),w=u({__name:"PdfViewer",setup(y){p(),c();const e=d();return l.openFileOS(e.pathAbsolute),(B,a)=>{const i=f("cv-button");return _(),h(b,null,[n(C,{pathArray:t(e).breadCrumbPathArray},{default:r(()=>[n(m,{icon:"icn-close",btnStyle:"soft",mini:"",onClick:t(e).exitViewer},null,8,["onClick"])]),_:1},8,["pathArray"]),F,k,v,n(i,{size:"default",onClick:a[0]||(a[0]=x=>t(l).openFileOS(t(e).pathAbsolute))},{default:r(()=>[s("Open PDF")]),_:1})],64)}}});export{w as default};

import{d as u,u as p,a as c,b as d,g as l,h as _,i as f,j as m,p as s,q as i,I as h,k as o,x as e,v as n,t as b,F as y}from"./index-JP_UVOjq.js";import{B as k}from"./BreadCrumbs-EbsX7ctM.js";const C=e("h3",null,"File type not supported.",-1),x={class:"error-msg"},B=e("p",null,[n("This file will open in your system's default application instead."),e("br"),n("If nothing happens, click the button below.")],-1),S=e("br",null,null,-1),V=u({__name:"UnknownViewer",setup(w){p(),c();const t=d();return l.openFileOS(t.pathAbsolute),(A,a)=>{const r=_("cv-button");return f(),m(y,null,[s(k,{pathArray:o(t).breadCrumbPathArray},{default:i(()=>[s(h,{icon:"icn-close",btnStyle:"soft",mini:"",onClick:o(t).exitViewer},null,8,["onClick"])]),_:1},8,["pathArray"]),C,e("p",x,[n(" We don't support "),e("i",null,"."+b(o(t).ext),1),n(" files. ")]),B,S,s(r,{size:"default",onClick:a[0]||(a[0]=F=>o(l).openFileOS(o(t).pathAbsolute))},{default:i(()=>[n("Open file")]),_:1})],64)}}});export{V as default};

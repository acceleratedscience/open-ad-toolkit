import{d as m,b as d,e as c,i as r,j as i,p as l,q as p,K as f,k as a,l as _,t as C,n as h,F as k}from"./index-M2pymOVI.js";import{B as x}from"./BreadCrumbs-BJsAHiin.js";import{T as y}from"./TextBox-Dt-kVu0X.js";const B={key:1},j=m({__name:"JsonViewer",props:{filePath:{},data:{}},setup(u){const t=d(),n=u,o=c(()=>{if(n.data)return n.data;const e=t.data;return e!=null&&e.cacheId&&(e!=null&&e.mols)?e.mols:e}),s=c(()=>o.value?JSON.stringify(o.value,null,4):"");return(e,b)=>(r(),i(k,null,[l(x,{pathArray:a(t).breadCrumbPathArray},{default:p(()=>[l(f,{icon:"icn-close",btnStyle:"soft",mini:"",onClick:a(t).exitViewer},null,8,["onClick"])]),_:1},8,["pathArray"]),s.value?(r(),_(y,{key:0,textData:s.value},null,8,["textData"])):a(t).errCode?(r(),i("div",B,C(a(t).errCode),1)):h("",!0)],64))}});export{j as _};

import{d as B,b as E,L as W,r as x,Q as N,e as C,G as T,w as y,i as r,j as o,p as _,q as V,K as D,k as h,F as f,N as g,O as I,t as b,n as F,E as M}from"./index-SnNRlG3y.js";import{B as L}from"./BreadCrumbs-BCeGk8Z9.js";const P={key:1},$=B({__name:"TextViewer--bkp",props:{filePath:{},data:{}},setup(k){const u=E(),w=W(),v=k,i=x(null),d=x(0),{screenWidth:A}=N(w),p=C(()=>v.data?v.data:u.data),S=C(()=>{if(!p.value)return[];let t=String(p.value).split(`
`);if(console.log(t),!d.value)return t;const n=new Array;return t.forEach(a=>{for(n.push([]);a.length>d.value;){let e=a.slice(0,d.value+1),s=d.value;if(e[e.length-1]==" ")e=e.slice(0,e.length-1);else{let l=[...e.matchAll(/\s/g)],c=l.length>0?l[l.length-1].index:-1;c&&c>=0&&(s=c,e=e.slice(0,s))}n[n.length-1].push(e.trim()),a=a.slice(s)}n[n.length-1].push(a.trim())}),console.log(n),n});function m(){if(!i.value)return 0;const t=i.value.cloneNode(!0);t.innerHTML="",t.style.width=i.value.offsetWidth+"px",t.style.position="absolute",t.style.top="-1000px",t.style.opacity="0";const n=document.createElement("div"),a=document.createElement("span");n.appendChild(a),t.appendChild(n),document.body.appendChild(t);const e=i.value.children[0];e&&e.getAttributeNames().forEach(c=>{n.setAttribute(c,e.getAttribute(c))});let s=0;const l=91;for(;a.offsetWidth+l<i.value.offsetWidth;)a.innerText+="m",s++;document.body.removeChild(t),d.value=s-1}return T(()=>m()),y(p,()=>m()),y(A,()=>m()),(t,n)=>(r(),o(f,null,[_(L,{pathArray:h(u).breadCrumbPathArray},{default:V(()=>[_(D,{icon:"icn-close",btnStyle:"soft",mini:"",onClick:h(u).exitViewer},null,8,["onClick"])]),_:1},8,["pathArray"]),p.value?(r(),o("pre",{key:0,ref_key:"$textContent",ref:i,class:"text-content"},[(r(!0),o(f,null,g(S.value,(a,e)=>(r(),o(f,{key:e},[(r(!0),o(f,null,g(a,(s,l)=>(r(),o("div",{key:l,class:I({"new-line":l===0})},b(s),3))),128))],64))),128))],512)):h(u).errCode?(r(),o("div",P,b(h(u).errCode),1)):F("",!0)],64))}}),G=M($,[["__scopeId","data-v-57591760"]]);export{G as default};

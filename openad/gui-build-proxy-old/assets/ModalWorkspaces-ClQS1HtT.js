import{d as S,L as W,M as x,r as i,G as y,g as a,h as u,i as b,l as h,q as n,v as c,p as C,k as M}from"./index-M2pymOVI.js";const g=S({__name:"ModalWorkspaces",emits:["mounted"],setup(V,{emit:d}){const m=W(),r=x(),p=d,l=i([" "]),s=i(" ");y(()=>{p("mounted"),v()});async function f(){if(!a)return;const{status:o,statusText:e}=await a.setWorkspace(s.value);o!==200?console.error(e):(m.setWorkspace(s.value),r.hide())}async function v(){const o=await k();if(o){const{all:e,active:t}=o;e&&(l.value=e),t&&(s.value=t)}else console.error("Failed to load workspaces")}async function k(){if(!a)return;const{status:o,data:e,statusText:t}=await a.getWorkspaces();if(o!==200)console.error(t);else return e}return(o,e)=>{const t=u("cv-dropdown"),_=u("cv-modal");return b(),h(_,{visible:M(r).visible,size:"xs",onPrimaryClick:f},{title:n(()=>[c("Switch workspace")]),content:n(()=>[C(t,{modelValue:s.value,"onUpdate:modelValue":e[0]||(e[0]=w=>s.value=w),items:l.value},null,8,["modelValue","items"])]),"secondary-button":n(()=>[c("Cancel")]),"primary-button":n(()=>[c("Switch")]),_:1},8,["visible"])}}});export{g as default};

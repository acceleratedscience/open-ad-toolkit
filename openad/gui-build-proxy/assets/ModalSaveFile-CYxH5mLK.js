import{d as R,L as h,R as z,r as b,c as x,H as A,g as T,h as s,k as r,p as t,q as u,v as c,t as v,j as w,n as p,i as f,F as y,A as C,m as D}from"./index-B1am-EuC.js";import W from"./FileBrowser-CcG7b_Pp.js";const $={class:"flex-wrap-v"},G={class:"input-pad"},J={class:"file-ext"},Y=R({__name:"ModalSaveFile",emits:["mounted"],setup(K,{emit:O}){const n=h(),B=z(),U=O,g=b(!1),E=[{val:"smol.json",disp:"smol.json"},{val:"sdf",disp:"SDF"},{val:"csv",disp:"CSV"},{val:"mol",disp:"MOL"},{val:"smi",disp:"SMILES"}],I=[{val:"mmol.json",disp:"mmol.json"},{val:"cif",disp:"CIF"}],L=[{val:"mmol.json",disp:"mmol.json"},{val:"pdb",disp:"PDB"},{val:"cif",disp:"CIF"}],P=[{val:"molset.json",disp:"molset.json"},{val:"sdf",disp:"SDF"},{val:"csv",disp:"CSV"},{val:"smi",disp:"SMILES"}],a=b(null),j=b(null),F=b(null),o=x(()=>n.data??{}),m=x({get:()=>o.value.filename??"",set:d=>{F.value=null,n.setData({filename:d})}}),_=x({get:()=>o.value.path??"",set:d=>n.setData({path:d})}),M=x(()=>a.value||"");A(()=>{U("mounted"),o.value.dataType=="smol"?a.value="smol.json":["cif","pdb"].includes(o.value.dataType)?a.value="mmol.json":o.value.dataType=="molset"&&(a.value="molset.json"),setTimeout(()=>{var d;j.value=document.querySelector("#modal-save-file input.bx--text-input"),(d=j.value)==null||d.select()},350)});async function q(){var V;if(!m.value){F.value="Filename required.",(V=j.value)==null||V.focus();return}const l=(_.value?_.value+"/":"")+m.value+"."+M.value;n.onSubmit&&n.onSubmit({destinationPath:l,ext:M.value,srcDataType:o.value.dataType}),g.value=!0}async function N(){!g.value&&n.onCancel&&n.onCancel()}return(d,l)=>{const V=T("cv-text-input"),S=T("cv-dropdown-item"),k=T("cv-dropdown"),H=T("cv-modal");return s(),r(H,{id:"modal-save-file",visible:w(n).visible,size:"lg",onPrimaryClick:q,onModalHidden:N},{title:t(()=>[u("Save to "),c("b",null,v(w(B).workspace),1),u(" Workspace")]),content:t(()=>[c("div",$,[c("div",G,[p(V,{modelValue:m.value,"onUpdate:modelValue":l[0]||(l[0]=e=>m.value=e),class:"ip-filename","hide-label":!0,"invalid-message":F.value},null,8,["modelValue","invalid-message"]),c("div",J,[c("span",null,v(m.value),1),u("."+v(M.value),1)]),o.value.exportOptions?(s(),f(y,{key:0},[o.value.dataType=="smol"?(s(),r(k,{key:0,modelValue:a.value,"onUpdate:modelValue":l[1]||(l[1]=e=>a.value=e),class:"dd-output-options"},{default:t(()=>[(s(),f(y,null,C(E,(e,i)=>p(S,{key:i,value:e.val},{default:t(()=>[u(v(e.disp),1)]),_:2},1032,["value"])),64))]),_:1},8,["modelValue"])):o.value.dataType=="cif"?(s(),r(k,{key:1,modelValue:a.value,"onUpdate:modelValue":l[2]||(l[2]=e=>a.value=e),class:"dd-output-options"},{default:t(()=>[(s(),f(y,null,C(I,(e,i)=>p(S,{key:i,value:e.val},{default:t(()=>[u(v(e.disp),1)]),_:2},1032,["value"])),64))]),_:1},8,["modelValue"])):o.value.dataType=="pdb"?(s(),r(k,{key:2,modelValue:a.value,"onUpdate:modelValue":l[3]||(l[3]=e=>a.value=e),class:"dd-output-options"},{default:t(()=>[(s(),f(y,null,C(L,(e,i)=>p(S,{key:i,value:e.val},{default:t(()=>[u(v(e.disp),1)]),_:2},1032,["value"])),64))]),_:1},8,["modelValue"])):o.value.dataType=="molset"?(s(),r(k,{key:3,modelValue:a.value,"onUpdate:modelValue":l[4]||(l[4]=e=>a.value=e),class:"dd-output-options"},{default:t(()=>[(s(),f(y,null,C(P,(e,i)=>p(S,{key:i,value:e.val},{default:t(()=>[u(v(e.disp),1)]),_:2},1032,["value"])),64))]),_:1},8,["modelValue"])):D("",!0)],64)):D("",!0)]),p(W,{isModal:!0,modelValue:_.value,"onUpdate:modelValue":l[5]||(l[5]=e=>_.value=e)},null,8,["modelValue"])])]),"secondary-button":t(()=>[u("Cancel")]),"primary-button":t(()=>[u("Save")]),_:1},8,["visible"])}}});export{Y as default};

Clazz.declarePackage ("J.adapter.readers.quantum");
Clazz.load (["J.adapter.readers.quantum.AdfReader"], "J.adapter.readers.quantum.AmsReader", null, function () {
c$ = Clazz.declareType (J.adapter.readers.quantum, "AmsReader", J.adapter.readers.quantum.AdfReader);
Clazz.overrideMethod (c$, "initializeReader", 
function () {
this.isADF = false;
});
});

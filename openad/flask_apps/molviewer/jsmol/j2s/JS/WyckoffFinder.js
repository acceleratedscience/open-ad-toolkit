Clazz.declarePackage ("JS");
Clazz.load (["java.util.Hashtable", "JU.V3"], "JS.WyckoffFinder", ["JU.JSJSONParser", "$.Measure", "$.P3", "$.P4", "$.PT", "$.Rdr", "JU.SimpleUnitCell", "JV.FileManager"], function () {
c$ = Clazz.decorateAsClass (function () {
this.positions = null;
this.centerings = null;
Clazz.instantialize (this, arguments);
}, JS, "WyckoffFinder");
Clazz.makeConstructor (c$, 
function () {
});
Clazz.defineMethod (c$, "getWyckoffFinder", 
function (vwr, sgname) {
var helper = JS.WyckoffFinder.helpers.get (sgname);
if (helper == null) {
var itno = JU.PT.parseInt (JU.PT.split (sgname, ":")[0]);
if (itno >= 1 && itno <= 230) {
var resource = this.getResource (vwr, "ita_" + itno + ".json");
if (resource != null) {
var its = resource.get ("its");
if (its != null) {
for (var i = its.size (); --i >= 0; ) {
var map = its.get (i);
if (sgname.equals (map.get ("itaFull"))) {
JS.WyckoffFinder.helpers.put (sgname, helper =  new JS.WyckoffFinder (map));
return helper;
}}
}}}}if (helper == null) {
if (JS.WyckoffFinder.nullHelper == null) JS.WyckoffFinder.nullHelper =  new JS.WyckoffFinder (null);
JS.WyckoffFinder.helpers.put (sgname, JS.WyckoffFinder.nullHelper);
}return helper;
}, "JV.Viewer,~S");
Clazz.makeConstructor (c$, 
 function (map) {
if (map != null) {
var wpos = map.get ("wpos");
this.positions = wpos.get ("pos");
var cent = wpos.get ("cent");
if (cent != null) {
this.centerings =  new Array (cent.size ());
for (var i = cent.size (); --i >= 0; ) {
this.centerings[i] = JS.WyckoffFinder.toPoint (cent.get (i));
}
}}}, "java.util.Map");
Clazz.defineMethod (c$, "getWyckoffPosition", 
function (p) {
if (this.positions == null) return "?";
for (var i = this.positions.size (); --i >= 0; ) {
var map = this.positions.get (i);
if (i == 0) {
return map.get ("label");
}var coords = map.get ("coord");
for (var c = 0, n = coords.size (); c < n; c++) {
if (JS.WyckoffFinder.getWyckoffCoord (coords, c).contains (p, this.centerings)) {
return map.get ("label");
}}
}
return "?";
}, "JU.P3");
Clazz.defineMethod (c$, "findPositionFor", 
function (p, letter) {
if (this.positions == null) return null;
for (var i = this.positions.size (); --i >= 0; ) {
var map = this.positions.get (i);
if (map.get ("label").equals (letter)) {
var coords = map.get ("coord");
if (coords != null) JS.WyckoffFinder.getWyckoffCoord (coords, 0).set (p);
return p;
}}
return null;
}, "JU.P3,~S");
c$.getWyckoffCoord = Clazz.defineMethod (c$, "getWyckoffCoord", 
 function (coords, c) {
var coord = coords.get (c);
if (Clazz.instanceOf (coord, String)) {
coords.set (c, coord =  new JS.WyckoffFinder.WyckoffPos (coord));
}return coord;
}, "JU.Lst,~N");
Clazz.defineMethod (c$, "getResource", 
 function (vwr, resource) {
try {
var r = JV.FileManager.getBufferedReaderForResource (vwr, this, "JS/", "sg/json/" + resource);
var data =  new Array (1);
if (JU.Rdr.readAllAsString (r, 2147483647, false, data, 0)) {
return  new JU.JSJSONParser ().parse (data[0], true);
}} catch (e) {
System.err.println (e.getMessage ());
}
return null;
}, "JV.Viewer,~S");
c$.toPoint = Clazz.defineMethod (c$, "toPoint", 
function (xyz) {
var s = JU.PT.split (xyz, ",");
return JU.P3.new3 (JU.PT.parseFloatFraction (s[0]), JU.PT.parseFloatFraction (s[1]), JU.PT.parseFloatFraction (s[2]));
}, "~S");
Clazz.pu$h(self.c$);
c$ = Clazz.decorateAsClass (function () {
this.point = null;
this.line = null;
this.plane = null;
this.type = 0;
this.xyz = null;
Clazz.instantialize (this, arguments);
}, JS.WyckoffFinder, "WyckoffPos");
c$.unitize = Clazz.defineMethod (c$, "unitize", 
 function (a) {
JU.SimpleUnitCell.unitizeDim (3, a);
return a;
}, "JU.P3");
Clazz.makeConstructor (c$, 
function (a) {
this.create (a);
}, "~S");
Clazz.defineMethod (c$, "create", 
 function (a) {
var b = JU.PT.split (a, ",");
var c = 0;
for (var d = 0; d < 3; d++) {
if (b[d].indexOf ('x') >= 0) {
c |= 1;
}if (b[d].indexOf ('y') >= 0) {
c |= 2;
}if (b[d].indexOf ('z') >= 0) {
c |= 4;
}}
var e;
var f;
var g;
switch (c) {
case 0:
this.type = 1;
this.point = JS.WyckoffFinder.toPoint (a);
break;
case 1:
case 2:
case 4:
this.type = 2;
e = JS.WyckoffFinder.WyckoffPos.ptFor (a, 0, 0, 0);
f = JS.WyckoffFinder.WyckoffPos.ptFor (a, 1, 1.27, 1.64);
f.sub2 (f, e);
f.normalize ();
this.point = JS.WyckoffFinder.WyckoffPos.unitize (JU.P3.newP (e));
this.line = JU.V3.newV (f);
break;
case 3:
case 5:
case 6:
this.type = 3;
e = JS.WyckoffFinder.WyckoffPos.ptFor (a, 0, 0, 0);
f = JS.WyckoffFinder.WyckoffPos.ptFor (a, 1.23, 1.47, 1.86);
g = JS.WyckoffFinder.WyckoffPos.ptFor (a, 0.1, 0.2, 0.3);
this.plane = JU.Measure.getPlaneThroughPoints (e, f, g, null, null,  new JU.P4 ());
break;
case 7:
break;
}
}, "~S");
c$.ptFor = Clazz.defineMethod (c$, "ptFor", 
 function (a, b, c, d) {
var e = JU.PT.split (a, ",");
var f = JS.WyckoffFinder.WyckoffPos.decodeXYZ (e[0], b, c, d);
var g = JS.WyckoffFinder.WyckoffPos.decodeXYZ (e[1], b, c, d);
var h = JS.WyckoffFinder.WyckoffPos.decodeXYZ (e[2], b, c, d);
return JU.P3.new3 (f, g, h);
}, "~S,~N,~N,~N");
c$.decodeXYZ = Clazz.defineMethod (c$, "decodeXYZ", 
 function (a, b, c, d) {
a = JU.PT.rep (a, "-", "+-");
a = JU.PT.rep (a, "x", "*x");
a = JU.PT.rep (a, "y", "*y");
a = JU.PT.rep (a, "z", "*z");
a = JU.PT.rep (a, "-*", "-");
a = JU.PT.rep (a, "+*", "+");
var e = 0;
var f = JU.PT.split (a, "+");
for (var g = f.length; --g >= 0; ) {
a = f[g];
if (a.length == 0) continue;
if (a.indexOf ('.') >= 0) {
e += JU.PT.parseFloat (a);
continue;
}var h = 0;
var i = 0;
var j = 0;
var k = 1;
switch (a.charAt (0)) {
case '-':
k = -1;
case '*':
j++;
break;
}
for (var l = a.length; --l >= j; ) {
var m = a.charAt (l);
switch (m) {
case 'x':
h = b;
break;
case 'y':
h = c;
break;
case 'z':
h = d;
break;
case '/':
h = 1 / h;
case '*':
k *= h;
h = 0;
break;
default:
var n = "0123456789".indexOf (m);
if (n < 0) System.err.println ("WH ????");
if (h == 0) {
h = n;
} else {
i = (i == 0 ? 10 : i * 10);
h += i * n;
}break;
}
}
e += k * h;
}
return e;
}, "~S,~N,~N,~N");
Clazz.defineMethod (c$, "contains", 
function (a, b) {
if (this.containsPt (a)) return true;
var c =  new JU.P3 ();
if (b != null) for (var d = b.length; --d >= 0; ) {
c.add2 (a, b[d]);
JS.WyckoffFinder.WyckoffPos.unitize (c);
if (this.containsPt (c)) return true;
}
return false;
}, "JU.P3,~A");
Clazz.defineMethod (c$, "containsPt", 
function (a) {
var b = 1;
switch (this.type) {
case 1:
b = a.distance (this.point);
break;
case 2:
var c = JU.P3.newP (a);
JU.Measure.projectOntoAxis (c, this.point, this.line, JS.WyckoffFinder.WyckoffPos.vtemp1);
b = a.distance (c);
break;
case 3:
b = JU.Measure.distanceToPlane (this.plane, a);
break;
}
return JS.WyckoffFinder.WyckoffPos.approx (b) == 0;
}, "JU.P3");
Clazz.defineMethod (c$, "set", 
function (a) {
switch (this.type) {
case 1:
a.setT (this.point);
break;
case 2:
JU.Measure.projectOntoAxis (a, this.point, this.line, JS.WyckoffFinder.WyckoffPos.vtemp1);
break;
case 3:
JU.Measure.getPlaneProjection (a, this.plane, JS.WyckoffFinder.WyckoffPos.vtemp1, JS.WyckoffFinder.WyckoffPos.vtemp1);
a.setT (JS.WyckoffFinder.WyckoffPos.vtemp1);
break;
}
}, "JU.P3");
c$.approx = Clazz.defineMethod (c$, "approx", 
 function (a) {
return JU.PT.approx (a, 1000);
}, "~N");
Clazz.defineStatics (c$,
"TYPE_POINT", 1,
"TYPE_LINE", 2,
"TYPE_PLANE", 3);
c$.vtemp1 = c$.prototype.vtemp1 =  new JU.V3 ();
c$ = Clazz.p0p ();
c$.helpers = c$.prototype.helpers =  new java.util.Hashtable ();
Clazz.defineStatics (c$,
"nullHelper", null);
});

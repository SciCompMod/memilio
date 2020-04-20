function fun1(x) {
  return Math.sin(x);
}
function fun2(x) {
  return Math.cos(3 * x);
}

function draw() {
  var canvas = document.getElementById('canvas');
  if (null == canvas || !canvas.getContext) return;

  var axes = {},
    ctx = canvas.getContext('2d');
  axes.x0 = 0.5 + 0.5 * canvas.width; // x0 pixels from left to x=0
  axes.y0 = 0.5 + 0.5 * canvas.height; // y0 pixels from top to y=0
  axes.scale = 40; // 40 pixels from x=0 to x=1
  axes.doNegativeX = true;

  showAxes(ctx, axes);
  funGraph(ctx, axes, fun1, 'rgb(11,153,11)', 1);
  funGraph(ctx, axes, fun2, 'rgb(66,44,255)', 2);
}

function funGraph(ctx, axes, func, color, thick) {
  var xx,
    yy,
    dx = 4,
    x0 = axes.x0,
    y0 = axes.y0,
    scale = axes.scale;
  var iMax = Math.round((ctx.canvas.width - x0) / dx);
  var iMin = axes.doNegativeX ? Math.round(-x0 / dx) : 0;
  ctx.beginPath();
  ctx.lineWidth = thick;
  ctx.strokeStyle = color;

  for (var i = iMin; i <= iMax; i++) {
    xx = dx * i;
    yy = scale * func(xx / scale);
    if (i == iMin) ctx.moveTo(x0 + xx, y0 - yy);
    else ctx.lineTo(x0 + xx, y0 - yy);
  }
  ctx.stroke();
}

function showAxes(ctx, axes) {
  var x0 = axes.x0,
    w = ctx.canvas.width;
  var y0 = axes.y0,
    h = ctx.canvas.height;
  var xmin = axes.doNegativeX ? 0 : x0;
  ctx.beginPath();
  ctx.strokeStyle = 'rgb(128,128,128)';
  ctx.moveTo(xmin, y0);
  ctx.lineTo(w, y0); // X axis
  ctx.moveTo(x0, 0);
  ctx.lineTo(x0, h); // Y axis
  ctx.stroke();
}

// Week number according to the ISO-8601 standard
Date.prototype.getWeek = function () {
  var target = new Date(this.valueOf());
  var dayNr = (this.getDay() + 6) % 7;
  target.setDate(target.getDate() - dayNr + 3);
  var firstThursday = target.valueOf();
  target.setMonth(0, 1);
  if (target.getDay() != 4) {
    target.setMonth(0, 1 + ((4 - target.getDay() + 7) % 7));
  }
  return 1 + Math.ceil((firstThursday - target) / 604800000);
};

function showFirstActionInterval(target) {
  var elem_div = document.getElementById('action' + target);
  var elem_input = document.getElementById('input_action' + target);

  //  + "

  if (elem_input.checked) {
    show('action' + target, 0);
  } else {
    hideAll('action' + target + '_');
  }
}

function show(target_name, number) {
  var elem = document.getElementById(target_name + '_' + number);

  if (elem.style.display == 'none') {
    elem.style.display = 'block';
  }
}

function hideAll(target_name) {
  for (var i = 0; i < 4; i++) {
    var elem = document.getElementById(target_name + '' + i);

    elem.style.display = 'none';
  }
}

var actions = new Array(
  'HomeOffice Empfehlung',
  'HomeOffice Anordnung',
  'Schlie&szlig;ung von Restaurant und Bars',
  'Schlie&szlig;ung von Gesch&auml;ften',
  'Schulschlie&szlig;ungen',
  'Kontaktverbot',
  'Ausgangssperre',
  'Einreiseverbote',
  'Systematisierung von Tests'
);

function init() {
  document
    .getElementById('button')
    .addEventListener('click', changeHTMLContent);
}

function changeHTMLContent() {
  var neu = 'neuer <b>fetter</b> Text';
  document.getElementById('absatz').innerHTML = neu;
}

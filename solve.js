// Newton (NOT Lagrange) polynomial reconstruction with exact BigInt rationals.
// Usage: node solve.js <input.json> [--coeffs]

const fs = require("fs");

// -------------------- BigInt Rational --------------------
class Rational {
  constructor(num, den = 1n) {
    if (den === 0n) throw new Error("Denominator cannot be zero");
    // normalize sign
    if (den < 0n) {
      num = -num;
      den = -den;
    }
    const g = gcd(abs(num), den);
    this.n = num / g;
    this.d = den / g;
  }
  static fromBigInt(x) { return new Rational(BigInt(x), 1n); }
  static zero() { return new Rational(0n, 1n); }
  static one() { return new Rational(1n, 1n); }

  add(other) { return new Rational(this.n * other.d + other.n * this.d, this.d * other.d); }
  sub(other) { return new Rational(this.n * other.d - other.n * this.d, this.d * other.d); }
  mul(other) { return new Rational(this.n * other.n, this.d * other.d); }
  div(other) {
    if (other.n === 0n) throw new Error("Division by zero");
    return new Rational(this.n * other.d, this.d * other.n);
  }
  neg() { return new Rational(-this.n, this.d); }
  isInteger() { return this.d === 1n; }
  toString() { return this.isInteger() ? this.n.toString() : `${this.n.toString()}/${this.d.toString()}`; }
}

function abs(x) { return x < 0n ? -x : x; }
function gcd(a, b) {
  while (b !== 0n) { const t = a % b; a = b; b = t; }
  return a < 0n ? -a : a;
}

// -------------------- Base parsing to BigInt --------------------
function parseInBaseToBigInt(str, base) {
  const b = BigInt(base);
  const s = str.trim().toLowerCase();
  let val = 0n;
  for (let i = 0; i < s.length; i++) {
    const ch = s[i];
    let digit;
    if (ch >= '0' && ch <= '9') digit = BigInt(ch.charCodeAt(0) - 48);
    else if (ch >= 'a' && ch <= 'z') digit = BigInt(ch.charCodeAt(0) - 87); // 'a'->10
    else throw new Error(`Invalid digit '${ch}'`);
    if (digit >= b) throw new Error(`Digit '${ch}' out of range for base ${base}`);
    val = val * b + digit;
  }
  return val;
}

// -------------------- Newton Divided Differences --------------------
function newtonDividedDifferences(xs, ys) {
  // xs: Rational[], ys: Rational[]
  const n = xs.length;
  const dd = Array.from({ length: n }, () => Array(n).fill(Rational.zero()));
  for (let i = 0; i < n; i++) dd[i][0] = ys[i];
  for (let j = 1; j < n; j++) {
    for (let i = 0; i < n - j; i++) {
      const numerator = dd[i + 1][j - 1].sub(dd[i][j - 1]);
      const denom = xs[i + j].sub(xs[i]);
      dd[i][j] = numerator.div(denom);
    }
  }
  // Return the top row c0..c_{n-1}
  return dd[0];
}

// Evaluate Newton form at x=0 directly (stable/no expansion)
function evaluateAtZeroNewton(xs, coeffs) {
  // coeffs: c0..c_{m}; xs: x0..x_{m-1}
  // f(0) = c0 + c1*(0 - x0) + c2*(0 - x0)(0 - x1)+ ...
  let result = coeffs[0];
  let prod = Rational.one();
  for (let j = 1; j < coeffs.length; j++) {
    prod = prod.mul(Rational.zero().sub(xs[j - 1])); // (0 - x_{j-1})
    result = result.add(coeffs[j].mul(prod));
  }
  return result;
}

// Optional: convert Newton form to standard polynomial coefficients a0..am
// Using: Start with P(x)=c_m, then for i=m-1..0: P(x)=c_i + (x - x_i)P(x)
function newtonToStandardCoeffs(xs, coeffs) {
  let poly = [coeffs[coeffs.length - 1]]; // degree 0
  for (let i = coeffs.length - 2; i >= 0; i--) {
    const xi = xs[i];
    const next = Array(poly.length + 1).fill(Rational.zero());
    // multiply by x: shift
    for (let j = 0; j < poly.length; j++) next[j + 1] = next[j + 1].add(poly[j]);
    // subtract xi * poly (to constant..)
    for (let j = 0; j < poly.length; j++) next[j] = next[j].sub(xi.mul(poly[j]));
    // add c_i to constant term
    next[0] = next[0].add(coeffs[i]);
    poly = next;
  }
  return poly; // a0, a1, ..., am
}

// -------------------- Main --------------------
function main() {
  const showCoeffs = process.argv.includes("--coeffs");
  const file = process.argv[2];
  if (!file) {
    console.error("Usage: node solve.js <input.json> [--coeffs]");
    process.exit(1);
  }

  const raw = fs.readFileSync(file, "utf8");
  const data = JSON.parse(raw);

  const k = Number(data.keys.k); // degree m = k-1
  // Build points (x, y) and convert y from its base to BigInt
  const points = [];
  for (const key of Object.keys(data)) {
    if (key === "keys") continue;
    const base = Number(data[key].base);
    const valStr = data[key].value;
    const y = parseInBaseToBigInt(valStr, base);
    const x = BigInt(key);
    points.push({ x, y });
  }

  // Sort by x, pick first k points
  points.sort((a, b) => (a.x < b.x ? -1 : a.x > b.x ? 1 : 0));
  const chosen = points.slice(0, k);

  // Build xs, ys as Rational
  const xs = chosen.map(p => new Rational(p.x, 1n));
  const ys = chosen.map(p => new Rational(p.y, 1n));

  // Newton coefficients (top row of divided difference table)
  const coeffs = newtonDividedDifferences(xs, ys);

  // Secret is f(0) (constant term)
  const secret = evaluateAtZeroNewton(xs, coeffs);

  console.log("Using points (x: y_dec):");
  chosen.forEach(p => console.log(`  ${p.x.toString()}: ${p.y.toString()}`));
  console.log(`k = ${k} -> degree m = ${k - 1}`);
  console.log(`Secret f(0): ${secret.toString()}`);

  if (showCoeffs) {
    const standard = newtonToStandardCoeffs(xs, coeffs);
    console.log("\nPolynomial coefficients (a0 + a1*x + ... + am*x^m):");
    standard.forEach((c, i) => {
      console.log(`  a${i} = ${c.toString()}`);
    });
  }
}

main();

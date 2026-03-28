/*
** benchmark_gauss.e
**
** dccelib Performance Benchmark — GAUSS
**
** Run from repo root:
**   C:\gauss26\tgauss.exe -b -nj validation\benchmark_gauss.e
*/

new;
library dccelib;

B_REPS  = 5;
OUTFILE = __FILE_DIR $+ "benchmark_gauss_results.csv";

// -----------------------------------------------------------------------
// Simulate a balanced panel (program-level, no local)
// -----------------------------------------------------------------------
proc (1) = _sim_panel(N, T_len, seed);
    local f, ids, yrs, x1, x2, fe, fac, y;
    rndseed seed;
    f    = rndn(T_len, 1);
    ids  = reshape(seqa(1,1,N) .* ones(1,T_len), N*T_len, 1);
    yrs  = reshape(ones(N,1) .* seqa(1,1,T_len)', N*T_len, 1);
    x1   = rndn(N*T_len, 1) + 0.5 .* reshape(ones(N,1) .* f', N*T_len, 1);
    x2   = rndn(N*T_len, 1) + 0.3 .* reshape(ones(N,1) .* f', N*T_len, 1);
    fe   = reshape(rndn(N,1) .* ones(1,T_len), N*T_len, 1);
    fac  = reshape(ones(N,1) .* f', N*T_len, 1);
    y    = 0.5.*x1 + 0.3.*x2 + fe + fac + 0.5.*rndn(N*T_len, 1);
    retp(asDF(ids~yrs~y~x1~x2, "id", "year", "y", "x1", "x2"));
endp;

// -----------------------------------------------------------------------
// Time one estimator B times, return median seconds
// -----------------------------------------------------------------------
proc (1) = _time_mg(data, ctl, B);
    local times, i, t0;
    struct mgOut mgO;
    times = zeros(B, 1);
    i = 1;
    do while i <= B;
        t0      = hsec();
        mgO     = mg(data, ctl);
        times[i]= (hsec()-t0)/100;
        i = i+1;
    endo;
    retp(quantile(times, 0.5));
endp;

proc (1) = _time_cce(data, ctl, B);
    local times, i, t0;
    struct mgOut cceO;
    times = zeros(B, 1);
    i = 1;
    do while i <= B;
        t0      = hsec();
        cceO    = cce_mg(data, ctl);
        times[i]= (hsec()-t0)/100;
        i = i+1;
    endo;
    retp(quantile(times, 0.5));
endp;

proc (1) = _time_dcce(data, ctl, B);
    local times, i, t0;
    struct mgOut dcceO;
    times = zeros(B, 1);
    i = 1;
    do while i <= B;
        t0      = hsec();
        dcceO   = dcce_mg(data, ctl);
        times[i]= (hsec()-t0)/100;
        i = i+1;
    endo;
    retp(quantile(times, 0.5));
endp;

proc (1) = _time_cips(data, B);
    local times, i, t0, cs, cv;
    times = zeros(B, 1);
    i = 1;
    do while i <= B;
        t0      = hsec();
        { cs, cv } = cips(data, 1, 1, 0);
        times[i]= (hsec()-t0)/100;
        i = i+1;
    endo;
    retp(quantile(times, 0.5));
endp;

proc (1) = _time_west(data, ctl, B);
    local times, i, t0;
    struct westerlundOut wO;
    times = zeros(B, 1);
    i = 1;
    do while i <= B;
        t0      = hsec();
        wO      = westerlundTest(data, 0, 1, 0);
        times[i]= (hsec()-t0)/100;
        i = i+1;
    endo;
    retp(quantile(times, 0.5));
endp;

// -----------------------------------------------------------------------
// Load Penn data
// -----------------------------------------------------------------------
fname    = __FILE_DIR $+ "../examples/penn_world.dta";
penn     = packr(loadd(fname, ". + date($year, '%Y')"));
penn     = order(penn, "id"$|"year");
penn_reg = penn[., "id" "year" "log_rgdpo" "log_ck" "log_ngd"];

N_penn = rows(unique(penn[., "id"]));
T_penn = rows(unique(penn[., "year"]));

struct mgControl ctl;
ctl = mgControlCreate();
ctl.report = 0;

struct mgControl ctl2;

// Rename penn columns to match sim naming for cips
penn_y = setcolnames(penn[., "id" "year" "log_rgdpo"], "id"$|"year"$|"y");

print "=============================================================";
print "dccelib Benchmark — GAUSS";
print "=============================================================";
print "";

// -----------------------------------------------------------------------
// Results accumulator: use fh to write directly to CSV
// -----------------------------------------------------------------------
fh = fopen(OUTFILE, "w");
fputs(fh, "dataset,N,T,estimator,tool,time_sec\n");

// Helper: write one row
// (called inline below via fputs)

// -----------------------------------------------------------------------
// Penn World Tables
// -----------------------------------------------------------------------
sprintf("Dataset: Penn World Tables  (N=%d, T~%d)", N_penn, T_penn);

t = _time_mg(penn_reg, ctl, B_REPS);
sprintf("  MG:         %.4f s", t);
fputs(fh, "Penn WTables," $+ ntos(N_penn,5) $+ "," $+ ntos(T_penn,5) $+ ",MG,GAUSS/dccelib," $+ ntos(t,6) $+ "\n");

ctl.x_csa = penn[., "log_hc"];
t = _time_cce(penn_reg, ctl, B_REPS);
sprintf("  CCE-MG:     %.4f s", t);
fputs(fh, "Penn WTables," $+ ntos(N_penn,5) $+ "," $+ ntos(T_penn,5) $+ ",CCE-MG,GAUSS/dccelib," $+ ntos(t,6) $+ "\n");

ctl.y_lags  = 1;
ctl.cr_lags = 3;
t = _time_dcce(penn_reg, ctl, B_REPS);
sprintf("  DCCE-MG:    %.4f s", t);
fputs(fh, "Penn WTables," $+ ntos(N_penn,5) $+ "," $+ ntos(T_penn,5) $+ ",DCCE-MG,GAUSS/dccelib," $+ ntos(t,6) $+ "\n");

ctl2 = mgControlCreate();
ctl2.report = 0;
t = _time_cips(penn_y, B_REPS);
sprintf("  CIPS:       %.4f s", t);
fputs(fh, "Penn WTables," $+ ntos(N_penn,5) $+ "," $+ ntos(T_penn,5) $+ ",CIPS,GAUSS/dccelib," $+ ntos(t,6) $+ "\n");

t = _time_west(penn_reg, ctl2, B_REPS);
sprintf("  Westerlund: %.4f s", t);
fputs(fh, "Penn WTables," $+ ntos(N_penn,5) $+ "," $+ ntos(T_penn,5) $+ ",Westerlund,GAUSS/dccelib," $+ ntos(t,6) $+ "\n");

print "";

// -----------------------------------------------------------------------
// Simulated panels
// -----------------------------------------------------------------------
struct mgControl ctl3;
struct mgControl ctl4;

sim_N = 100 | 200;
sim_T = 50  | 100;
sim_names = "Simulated M" $| "Simulated L";

s = 1;
do while s <= 2;
    sN    = sim_N[s];
    sT    = sim_T[s];
    sname = sim_names[s];
    sdata = _sim_panel(sN, sT, 42+s);

    sprintf("Dataset: %s  (N=%d, T=%d)", sname, sN, sT);

    ctl3 = mgControlCreate();
    ctl3.report = 0;

    t = _time_mg(sdata, ctl3, B_REPS);
    sprintf("  MG:         %.4f s", t);
    fputs(fh, sname $+ "," $+ ntos(sN,5) $+ "," $+ ntos(sT,5) $+ ",MG,GAUSS/dccelib," $+ ntos(t,6) $+ "\n");

    t = _time_cce(sdata, ctl3, B_REPS);
    sprintf("  CCE-MG:     %.4f s", t);
    fputs(fh, sname $+ "," $+ ntos(sN,5) $+ "," $+ ntos(sT,5) $+ ",CCE-MG,GAUSS/dccelib," $+ ntos(t,6) $+ "\n");

    ctl3.y_lags  = 1;
    ctl3.cr_lags = 3;
    t = _time_dcce(sdata, ctl3, B_REPS);
    sprintf("  DCCE-MG:    %.4f s", t);
    fputs(fh, sname $+ "," $+ ntos(sN,5) $+ "," $+ ntos(sT,5) $+ ",DCCE-MG,GAUSS/dccelib," $+ ntos(t,6) $+ "\n");

    sdata_y = sdata[., "id" "year" "y"];
    t = _time_cips(sdata_y, B_REPS);
    sprintf("  CIPS:       %.4f s", t);
    fputs(fh, sname $+ "," $+ ntos(sN,5) $+ "," $+ ntos(sT,5) $+ ",CIPS,GAUSS/dccelib," $+ ntos(t,6) $+ "\n");

    ctl4 = mgControlCreate();
    ctl4.report = 0;
    t = _time_west(sdata, ctl4, B_REPS);
    sprintf("  Westerlund: %.4f s", t);
    fputs(fh, sname $+ "," $+ ntos(sN,5) $+ "," $+ ntos(sT,5) $+ ",Westerlund,GAUSS/dccelib," $+ ntos(t,6) $+ "\n");

    print "";
    s = s+1;
endo;

fh = close(fh);

print "Results saved to: " $+ OUTFILE;

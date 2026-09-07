// tblastn-lite - find the gene for a protein: 6-frame translated nucleotide search.
// Neighbourhood word seeds + ungapped X-drop extension + Karlin-Altschul E-values.
// Single file, no dependencies.  MIT.
#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

// ---------------------------------------------------------------- scoring ---
static const char* AA = "ARNDCQEGHILKMFPSTWYVBZX*";  // BLOSUM62 column order
static const signed char BLOSUM62[24][24] = {
{ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4},
{-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4},
{-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4},
{-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1,-4},
{ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4},
{-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4},
{-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4},
{ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4},
{-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4},
{-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4},
{-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4},
{-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4},
{-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4},
{-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1,-4},
{-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2,-4},
{ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4},
{ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0,-4},
{-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2,-4},
{-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1,-4},
{ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4},
{-2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4},
{-1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4},
{ 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4},
{-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1}};
static const int STOP = 23, XAA = 22;
// ungapped BLOSUM62 Karlin-Altschul parameters (NCBI values)
static const double KA_LAMBDA = 0.318, KA_K = 0.134;

static signed char aa2i[256];
static signed char nt2i[256];  // T,C,A,G -> 0..3 ; anything else -> -1
static const char* CODONS =
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
static signed char codon2aa[64];

static void init_tables() {
  memset(aa2i, XAA, sizeof aa2i);
  for (int i = 0; i < 24; i++) {
    aa2i[(unsigned char)AA[i]] = (signed char)i;
    aa2i[(unsigned char)tolower(AA[i])] = (signed char)i;
  }
  aa2i[(unsigned char)'U'] = aa2i[(unsigned char)'u'] = 4;   // selenocysteine -> C
  aa2i[(unsigned char)'O'] = aa2i[(unsigned char)'o'] = 11;  // pyrrolysine   -> K
  memset(nt2i, -1, sizeof nt2i);
  const char* b = "TCAG";
  for (int i = 0; i < 4; i++) {
    nt2i[(unsigned char)b[i]] = (signed char)i;
    nt2i[(unsigned char)tolower(b[i])] = (signed char)i;
  }
  nt2i[(unsigned char)'U'] = nt2i[(unsigned char)'u'] = 0;
  for (int i = 0; i < 64; i++) codon2aa[i] = aa2i[(unsigned char)CODONS[i]];
}

// ------------------------------------------------------------------- fasta ---
struct Rec { std::string id, desc, seq; };

static std::vector<Rec> read_fasta(const std::string& path) {
  std::ifstream in(path);
  if (!in) { std::cerr << "error: cannot open " << path << "\n"; exit(1); }
  std::vector<Rec> recs;
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    if (line[0] == '>') {
      recs.push_back(Rec());
      size_t sp = line.find_first_of(" \t");
      recs.back().id = line.substr(1, sp == std::string::npos ? sp : sp - 1);
      if (sp != std::string::npos) recs.back().desc = line.substr(sp + 1);
    } else if (!recs.empty()) {
      recs.back().seq += line;
    }
  }
  if (recs.empty()) { std::cerr << "error: no sequences in " << path << "\n"; exit(1); }
  return recs;
}

static std::string revcomp(const std::string& s) {
  static const char* FWD = "ACGTacgt";
  static const char* REV = "TGCAtgca";
  std::string r(s.rbegin(), s.rend());
  for (size_t i = 0; i < r.size(); i++) {
    const char* p = strchr(FWD, r[i]);
    r[i] = (p && *p) ? REV[p - FWD] : 'N';
  }
  return r;
}

// translate s[offset..] into amino-acid codes
static void translate(const std::string& s, int offset, std::vector<signed char>& out) {
  size_t n = s.size() > (size_t)offset ? (s.size() - offset) / 3 : 0;
  out.resize(n);
  for (size_t i = 0; i < n; i++) {
    int a = nt2i[(unsigned char)s[offset + 3 * i]];
    int b = nt2i[(unsigned char)s[offset + 3 * i + 1]];
    int c = nt2i[(unsigned char)s[offset + 3 * i + 2]];
    out[i] = (a < 0 || b < 0 || c < 0) ? (signed char)XAA : codon2aa[16 * a + 4 * b + c];
  }
}

// ------------------------------------------------------------------ search ---
struct Opts {
  double evalue = 10.0;
  int thresh = 13, xdrop = 16, threads = 0, maxhits = 500;
  bool show_dna = false, show_aln = false;
};

struct Hit {
  int qs, qe;        // 0-based query aa range, inclusive
  long long ss, se;  // 1-based nucleotide coords on the plus strand (se < ss if minus)
  int frame, score, ident, len, seqidx;
  double evalue, bits;
};

// neighbourhood 3-mer index: word -> query positions
struct QueryIndex { std::vector<std::vector<int> > hits; };

static QueryIndex build_index(const std::vector<signed char>& q, const Opts& o) {
  QueryIndex ix;
  ix.hits.assign(8000, std::vector<int>());
  int rowmax[24];
  for (int a = 0; a < 24; a++) {
    int m = -100;
    for (int b = 0; b < 20; b++) m = std::max(m, (int)BLOSUM62[a][b]);
    rowmax[a] = m;
  }
  for (size_t i = 0; i + 3 <= q.size(); i++) {
    int q0 = q[i], q1 = q[i + 1], q2 = q[i + 2];
    if (q0 == STOP || q1 == STOP || q2 == STOP) continue;
    for (int a = 0; a < 20; a++) {
      int s1 = BLOSUM62[q0][a];
      if (s1 + rowmax[q1] + rowmax[q2] < o.thresh) continue;
      for (int b = 0; b < 20; b++) {
        int s2 = s1 + BLOSUM62[q1][b];
        if (s2 + rowmax[q2] < o.thresh) continue;
        for (int c = 0; c < 20; c++)
          if (s2 + BLOSUM62[q2][c] >= o.thresh)
            ix.hits[(a * 20 + b) * 20 + c].push_back((int)i);
      }
    }
  }
  return ix;
}

struct Scan {
  const std::vector<signed char>* q;
  const QueryIndex* ix;
  const Opts* o;
  double search_space;
};

static void scan_frame(const Scan& S, const std::vector<signed char>& t, int frame,
                       int seqidx, long long nlen, std::vector<Hit>& out) {
  const std::vector<signed char>& q = *S.q;
  const Opts& o = *S.o;
  int qn = (int)q.size(), tn = (int)t.size();
  if (tn < 3) return;
  std::vector<int> diag_end((size_t)qn + tn + 1, -1);  // last subject pos covered, per diagonal
  for (int sp = 0; sp + 3 <= tn; sp++) {
    int c0 = t[sp], c1 = t[sp + 1], c2 = t[sp + 2];
    if (c0 >= 20 || c1 >= 20 || c2 >= 20) continue;
    const std::vector<int>& bucket = S.ix->hits[(c0 * 20 + c1) * 20 + c2];
    for (size_t bi = 0; bi < bucket.size(); bi++) {
      int qp = bucket[bi], d = sp - qp + qn;
      if (diag_end[d] >= sp) continue;  // already inside an HSP on this diagonal
      int s = 0, bestR = 0, offR = -1;
      for (int k = 0; qp + k < qn && sp + k < tn; k++) {
        s += BLOSUM62[q[qp + k]][t[sp + k]];
        if (s > bestR) { bestR = s; offR = k; }
        if (bestR - s > o.xdrop) break;
      }
      if (offR < 0) continue;
      s = 0;
      int bestL = 0, offL = 0;
      for (int k = 1; qp - k >= 0 && sp - k >= 0; k++) {
        s += BLOSUM62[q[qp - k]][t[sp - k]];
        if (s > bestL) { bestL = s; offL = k; }
        if (bestL - s > o.xdrop) break;
      }
      int score = bestL + bestR, qs = qp - offL, ss = sp - offL, len = offL + offR + 1;
      diag_end[d] = ss + len - 1;
      double bits = (KA_LAMBDA * score - log(KA_K)) / log(2.0);
      double ev = S.search_space * pow(2.0, -bits);
      if (ev > o.evalue) continue;
      int ident = 0;
      for (int k = 0; k < len; k++) ident += (q[qs + k] == t[ss + k] && q[qs + k] < 20);
      Hit h;
      h.qs = qs; h.qe = qs + len - 1; h.frame = frame; h.score = score;
      h.ident = ident; h.len = len; h.seqidx = seqidx; h.evalue = ev; h.bits = bits;
      if (frame > 0) {
        h.ss = (long long)(frame - 1) + 3LL * ss + 1;
        h.se = h.ss + 3LL * len - 1;
      } else {
        h.ss = nlen - ((long long)(-frame - 1) + 3LL * ss);
        h.se = h.ss - 3LL * len + 1;
      }
      out.push_back(h);
    }
  }
}

// -------------------------------------------------------------------- main ---
static std::string aa_string(const std::vector<signed char>& v, int from, int len) {
  std::string s(len, 'X');
  for (int i = 0; i < len; i++) s[i] = AA[v[from + i]];
  return s;
}

// aa offset within its own frame for a hit's subject start
static int hit_frame_offset(const Hit& h, long long nlen) {
  int off = abs(h.frame) - 1;
  return (int)((h.frame > 0 ? h.ss - 1 - off : nlen - h.ss - off) / 3);
}

static int selftest();

// double-clicked on Windows: no arguments and the console dies with the process,
// so ask for the two paths instead of flashing the usage text.
static bool g_pause = false;
static void pause_exit() {
  if (!g_pause) return;
  std::cerr << "\nPress Enter to close...";
  std::cin.clear();
  std::string s;
  std::getline(std::cin, s);
}

static std::string ask(const char* msg) {
  std::cerr << msg;
  std::string s;
  if (!std::getline(std::cin, s)) return "";
  while (!s.empty() && (s.back() == ' ' || s.back() == '\r' || s.back() == '"')) s.erase(s.size() - 1);
  while (!s.empty() && (s[0] == ' ' || s[0] == '"')) s.erase(0, 1);
  return s;  // drag-and-dropped paths arrive quoted
}

int main(int argc, char** argv) {
  init_tables();
  std::string qpath, dpath;
  Opts o;

  if (argc == 1) {
    g_pause = true;
    atexit(pause_exit);
    std::cerr << "tblastn-lite - find the DNA encoding a protein\n"
                 "(drag a file onto this window to paste its path)\n\n";
    qpath = ask("Protein FASTA : ");
    dpath = ask("Genome FASTA  : ");
    std::string e = ask("E-value cutoff [1e-5] : ");
    if (!e.empty()) o.evalue = atof(e.c_str());
    else o.evalue = 1e-5;
    std::string out = ask("Save results to [blank = this window] : ");
    o.show_aln = o.show_dna = true;
    if (!out.empty() && !freopen(out.c_str(), "w", stdout)) {
      std::cerr << "error: cannot write " << out << "\n";
      return 1;
    }
    std::cerr << "\nsearching...\n";
  }


  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    bool needs = (a == "-q" || a == "--query" || a == "-d" || a == "--db" || a == "-e" ||
                  a == "--evalue" || a == "-T" || a == "--thresh" || a == "-X" ||
                  a == "--xdrop" || a == "-t" || a == "--threads" || a == "-m" || a == "--max");
    std::string v;
    if (needs) {
      if (i + 1 >= argc) { std::cerr << "error: " << a << " needs a value\n"; return 1; }
      v = argv[++i];
    }
    if (a == "-q" || a == "--query") qpath = v;
    else if (a == "-d" || a == "--db") dpath = v;
    else if (a == "-e" || a == "--evalue") o.evalue = atof(v.c_str());
    else if (a == "-T" || a == "--thresh") o.thresh = atoi(v.c_str());
    else if (a == "-X" || a == "--xdrop") o.xdrop = atoi(v.c_str());
    else if (a == "-t" || a == "--threads") o.threads = atoi(v.c_str());
    else if (a == "-m" || a == "--max") o.maxhits = atoi(v.c_str());
    else if (a == "--dna") o.show_dna = true;
    else if (a == "--aln") o.show_aln = true;
    else if (a == "--selftest") return selftest();
    else if (a == "-h" || a == "--help") { qpath.clear(); break; }
    else { std::cerr << "error: unknown option " << a << "\n"; return 1; }
  }
  if (qpath.empty() || dpath.empty()) {
    std::cerr <<
      "tblastn-lite - find the DNA encoding a protein (6-frame translated search)\n\n"
      "usage: tblastn_lite -q protein.fasta -d genome.fasta [options]\n"
      "  -e FLOAT   E-value cutoff (default 10)\n"
      "  -T INT     neighbourhood word threshold (default 13; lower = more sensitive)\n"
      "  -X INT     ungapped X-drop (default 16)\n"
      "  -t INT     threads (default: all cores)\n"
      "  -m INT     max HSPs reported per query (default 500)\n"
      "  --aln      print the alignment\n"
      "  --dna      print the matching DNA as FASTA\n"
      "  --selftest run built-in checks\n\n"
      "columns: qid sid %id len mism frame qstart qend sstart send evalue bits\n";
    return 1;
  }

  std::vector<Rec> db = read_fasta(dpath), queries = read_fasta(qpath);
  long long dblen = 0;
  for (size_t i = 0; i < db.size(); i++) dblen += (long long)db[i].seq.size();
  std::vector<std::string> rc(db.size());
  for (size_t i = 0; i < db.size(); i++) rc[i] = revcomp(db[i].seq);

  int nthr = o.threads > 0 ? o.threads : (int)std::max(1u, std::thread::hardware_concurrency());
  printf("# tblastn-lite  db: %s (%d seqs, %lld bp)\n", dpath.c_str(), (int)db.size(), dblen);

  struct Job { int seq, frame; };
  std::vector<Job> jobs;
  for (size_t i = 0; i < db.size(); i++)
    for (int f = 1; f <= 3; f++) {
      Job a = {(int)i, f}, b = {(int)i, -f};
      jobs.push_back(a); jobs.push_back(b);
    }

  for (size_t qi = 0; qi < queries.size(); qi++) {
    const Rec& Q = queries[qi];
    std::vector<signed char> q(Q.seq.size());
    for (size_t i = 0; i < Q.seq.size(); i++) q[i] = aa2i[(unsigned char)Q.seq[i]];
    if (q.size() < 3) { std::cerr << "skipping short query " << Q.id << "\n"; continue; }
    QueryIndex ix = build_index(q, o);
    Scan S = {&q, &ix, &o, (double)q.size() * (double)dblen * 2.0};

    std::vector<Hit> hits;
    std::mutex mtx;
    std::atomic<size_t> next(0);
    std::vector<std::thread> pool;
    for (int i = 0; i < nthr; i++) pool.push_back(std::thread([&]() {
      std::vector<signed char> t;
      std::vector<Hit> local;
      for (size_t j = next++; j < jobs.size(); j = next++) {
        const Job& job = jobs[j];
        const std::string& s = job.frame > 0 ? db[job.seq].seq : rc[job.seq];
        translate(s, abs(job.frame) - 1, t);
        scan_frame(S, t, job.frame, job.seq, (long long)db[job.seq].seq.size(), local);
      }
      std::lock_guard<std::mutex> g(mtx);
      hits.insert(hits.end(), local.begin(), local.end());
    }));
    for (size_t i = 0; i < pool.size(); i++) pool[i].join();

    std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) {
      return a.evalue != b.evalue ? a.evalue < b.evalue : a.score > b.score;
    });
    if ((int)hits.size() > o.maxhits) hits.resize(o.maxhits);

    printf("# query: %s (%d aa)  hits: %d\n", Q.id.c_str(), (int)Q.seq.size(), (int)hits.size());
    printf("# qid\tsid\t%%id\tlen\tmism\tframe\tqstart\tqend\tsstart\tsend\tevalue\tbits\n");
    std::vector<signed char> t;
    for (size_t hi = 0; hi < hits.size(); hi++) {
      const Hit& h = hits[hi];
      printf("%s\t%s\t%.1f\t%d\t%d\t%+d\t%d\t%d\t%lld\t%lld\t%.2g\t%.1f\n",
             Q.id.c_str(), db[h.seqidx].id.c_str(), 100.0 * h.ident / h.len, h.len,
             h.len - h.ident, h.frame, h.qs + 1, h.qe + 1, h.ss, h.se, h.evalue, h.bits);
      if (!o.show_aln && !o.show_dna) continue;
      const std::string& s = h.frame > 0 ? db[h.seqidx].seq : rc[h.seqidx];
      int off = abs(h.frame) - 1;
      translate(s, off, t);
      int sa0 = hit_frame_offset(h, (long long)db[h.seqidx].seq.size());
      if (o.show_aln) {
        std::string qa = aa_string(q, h.qs, h.len), sa = aa_string(t, sa0, h.len);
        std::string mid(h.len, ' ');
        for (int i = 0; i < h.len; i++)
          mid[i] = qa[i] == sa[i] ? qa[i] : (BLOSUM62[q[h.qs + i]][t[sa0 + i]] > 0 ? '+' : ' ');
        const int W = 60;
        for (int i = 0; i < h.len; i += W) {
          int n = std::min(W, h.len - i);
          long long s0 = h.frame > 0 ? h.ss + 3LL * i : h.ss - 3LL * i;
          long long s1 = h.frame > 0 ? s0 + 3LL * n - 1 : s0 - 3LL * n + 1;
          printf("Query %7d  %.*s  %d\n              %.*s\nSbjct %7lld  %.*s  %lld\n\n",
                 h.qs + i + 1, n, qa.c_str() + i, h.qs + i + n, n, mid.c_str() + i, s0, n,
                 sa.c_str() + i, s1);
        }
      }
      if (o.show_dna)
        printf(">%s:%lld-%lld frame%+d\n%s\n", db[h.seqidx].id.c_str(), h.ss, h.se, h.frame,
               s.substr(off + 3 * sa0, 3 * h.len).c_str());
    }
    if (g_pause) std::cerr << Q.id << ": " << hits.size() << " hits\n";
  }
  fflush(stdout);
  return 0;
}

// ---------------------------------------------------------------- selftest ---
static int selftest() {
  init_tables();
  assert(revcomp("ACGTNacgt") == "acgtNACGT");
  std::vector<signed char> t;
  translate("ATGGCTTAA", 0, t);
  assert(t.size() == 3 && AA[t[0]] == 'M' && AA[t[1]] == 'A' && AA[t[2]] == '*');
  translate("GATGGCTTAA", 1, t);
  assert(AA[t[0]] == 'M');

  // plant a known protein in random DNA, then find it from both strands
  const char* prot = "MKAIVLTYNQEWFRSGHDPLLNQVKACGEMTIYRDWVSAHFPNGKLMDEQRTYIVAKLNG";
  std::string dna;
  for (const char* p = prot; *p; p++)
    for (int c = 0; c < 64; c++)
      if (CODONS[c] == *p) {
        dna += "TCAG"[c / 16]; dna += "TCAG"[(c / 4) % 4]; dna += "TCAG"[c % 4];
        break;
      }
  assert(dna.size() == 3 * strlen(prot));
  std::string pad;
  srand(1);
  for (int i = 0; i < 500; i++) pad += "ACGT"[rand() % 4];
  std::string genome = pad + dna + pad;

  std::vector<signed char> q(strlen(prot));
  for (size_t i = 0; i < q.size(); i++) q[i] = aa2i[(unsigned char)prot[i]];
  Opts o;
  QueryIndex ix = build_index(q, o);
  Scan S = {&q, &ix, &o, (double)q.size() * (double)genome.size() * 2.0};
  for (int strand = 0; strand < 2; strand++) {
    std::string g = strand ? revcomp(genome) : genome, r = revcomp(g);
    std::vector<Hit> hits;
    for (int f = 1; f <= 3; f++) {
      translate(g, f - 1, t);
      scan_frame(S, t, f, 0, (long long)g.size(), hits);
      translate(r, f - 1, t);
      scan_frame(S, t, -f, 0, (long long)g.size(), hits);
    }
    assert(!hits.empty());
    std::sort(hits.begin(), hits.end(),
              [](const Hit& a, const Hit& b) { return a.bits > b.bits; });
    const Hit& h = hits[0];
    assert(h.len == (int)q.size() && h.ident == (int)q.size());  // full-length exact hit
    assert(h.qs == 0 && h.qe == (int)q.size() - 1);
    assert(h.evalue < 1e-20);
    long long lo = std::min(h.ss, h.se), hi = std::max(h.ss, h.se);
    assert(hi - lo + 1 == (long long)dna.size());
    std::string got = g.substr(lo - 1, dna.size());
    assert(got == (h.frame > 0 ? dna : revcomp(dna)));
    assert(hit_frame_offset(h, (long long)g.size()) * 3 + abs(h.frame) - 1 ==
           (int)(h.frame > 0 ? lo - 1 : (long long)g.size() - hi));
  }
  printf("selftest OK\n");
  return 0;
}

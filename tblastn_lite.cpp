// tblastn-lite - command-line front end.  See tblastn_core.h for the search itself.
#include "tblastn_core.h"

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
  std::string qpath, dpath, acc;
  Opts o;

  if (argc == 1) {
    g_pause = true;
    atexit(pause_exit);
    std::cerr << "tblastn-lite - find the DNA encoding a protein\n"
                 "(drag a file onto this window to paste its path)\n\n";
    qpath = ask("Protein FASTA or proteome TSV : ");
    acc = ask("Accession (blank = all entries) : ");
    dpath = ask("Genome FASTA  : ");
    std::string e = ask("E-value cutoff [1e-5] : ");
    o.evalue = e.empty() ? 1e-5 : atof(e.c_str());
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
                  a == "--xdrop" || a == "-t" || a == "--threads" || a == "-m" ||
                  a == "--max" || a == "-a" || a == "--acc");
    std::string v;
    if (needs) {
      if (i + 1 >= argc) { std::cerr << "error: " << a << " needs a value\n"; return 1; }
      v = argv[++i];
    }
    if (a == "-q" || a == "--query") qpath = v;
    else if (a == "-a" || a == "--acc") acc = v;
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
      "       tblastn_lite -q proteome.tsv -a P00561 -d genome.fasta\n\n"
      "  -q FILE    protein FASTA, or a UniProt proteome TSV (Entry+Sequence columns)\n"
      "  -a ACC     one entry by accession or gene name (default: every entry)\n"
      "  -e FLOAT   E-value cutoff (default 10)\n"
      "  -T INT     neighbourhood word threshold (default 13; lower = more sensitive)\n"
      "  -X INT     ungapped X-drop (default 16)\n"
      "  -t INT     threads (default: all cores)\n"
      "  -m INT     max HSPs reported per query (default 500)\n"
      "  --aln      print the alignment\n"
      "  --dna      print the matching DNA as FASTA\n"
      "  --selftest run built-in checks\n\n"
      "run with no arguments for interactive prompts; columns are\n"
      "qid sid %id len mism frame qstart qend sstart send evalue bits\n";
    return 1;
  }

  Db db;
  std::vector<Rec> queries;
  try {
    db = load_db(dpath);
    queries = read_queries(qpath);
    if (!acc.empty()) {
      int i = find_accession(queries, acc);
      if (i < 0) throw std::runtime_error("no entry '" + acc + "' in " + qpath);
      queries = std::vector<Rec>(1, queries[i]);
    }
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
  printf("# tblastn-lite  db: %s (%d seqs, %lld bp)\n", dpath.c_str(), (int)db.recs.size(), db.len);
  for (size_t i = 0; i < queries.size(); i++) {
    std::vector<signed char> q = encode_protein(queries[i].seq);
    if (q.size() < 3) { std::cerr << "skipping short query " << queries[i].id << "\n"; continue; }
    if (!queries[i].desc.empty())
      printf("# %s  %s\n", queries[i].id.c_str(), queries[i].desc.c_str());
    std::vector<Hit> hits = search(q, db, o, o.threads);
    fputs(format_hits(queries[i].id, q, db, hits, o).c_str(), stdout);
    if (g_pause) std::cerr << queries[i].id << ": " << hits.size() << " hits\n";
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

  std::vector<signed char> q = encode_protein(prot);
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
  // the threaded driver must find the same thing as the raw frame scan
  Db db;
  db.recs.push_back(Rec());
  db.recs[0].id = "plant";
  db.recs[0].seq = genome;
  db.rc.push_back(revcomp(genome));
  db.len = (long long)genome.size();
  std::vector<Hit> hits = search(q, db, o, 4);
  assert(!hits.empty() && hits[0].ident == (int)q.size());
  assert(format_hits("p", q, db, hits, o).find("plant") != std::string::npos);

  // proteome TSV: columns found by name, not position, and accession/gene lookup
  const char* tmp = "tblastn_selftest.tsv";
  {
    std::ofstream f(tmp);
    f << "Length\tEntry\tJunk\tProtein names\tGene Names\tSequence\n"
      << "3\tP12345\t-\tFake protein\tfakA b0001\tMKA\n"
      << "4\tQ99999\t-\tOther\totherB\tMKAI\n"
      << "0\tR00000\t-\tNo sequence\tnoseq\t\n";
  }
  std::vector<Rec> tsv = read_queries(tmp);
  remove(tmp);
  assert(tsv.size() == 2);                       // the row with no sequence is dropped
  assert(tsv[0].id == "P12345" && tsv[0].seq == "MKA");
  assert(tsv[0].desc.substr(0, 4) == "fakA");    // first gene name, then protein name
  assert(find_accession(tsv, "p12345") == 0);    // accession, case-insensitive
  assert(find_accession(tsv, "otherB") == 1);    // or gene name
  assert(find_accession(tsv, "nope") == -1);
  printf("selftest OK\n");
  return 0;
}

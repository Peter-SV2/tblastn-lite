// tblastn-lite - plain Win32 GUI.  No toolkit, no runtime: pick a protein, pick a genome, Search.
// Build: g++ -O3 -std=c++17 -static -mwindows -o tblastn_gui tblastn_gui.cpp -lcomdlg32 -lshell32
#include <windows.h>
#include <commdlg.h>
#include <shellapi.h>

#include "tblastn_core.h"

enum {
  IDC_QPATH = 101, IDC_QPICK, IDC_ACC, IDC_DPATH, IDC_DPICK, IDC_EVAL,
  IDC_ALN, IDC_DNA, IDC_RUN, IDC_SAVE, IDC_OUT, IDC_STATUS,
  IDC_L1 = 201, IDC_L2, IDC_L3, IDC_L4
};
static const UINT WM_DONE = WM_APP + 1;  // wParam = hit count, lParam = new std::string*

static HWND hQ, hPick, hAcc, hD, hDPick, hE, hAln, hDna, hRun, hSave, hOut, hStatus;
static HWND hL1, hL2, hL3, hL4;
static HFONT hMono, hUi;
static std::vector<Rec> g_queries;   // entries of the loaded FASTA / proteome TSV
static std::string g_loaded;         // path g_queries came from
static std::string g_text;           // full result, unabridged, for Save
static bool g_busy = false;

static std::string get_text(HWND h) {
  int n = GetWindowTextLengthA(h);
  std::string s((size_t)n, '\0');
  if (n) GetWindowTextA(h, &s[0], n + 1);
  return s;
}

// EDIT controls want CRLF
static std::string to_crlf(const std::string& s) {
  std::string o;
  o.reserve(s.size() + s.size() / 20);
  for (size_t i = 0; i < s.size(); i++) {
    if (s[i] == '\n') o += '\r';
    o += s[i];
  }
  return o;
}

static std::string pick_file(HWND owner, const char* filter, bool save) {
  char buf[MAX_PATH] = "";
  OPENFILENAMEA ofn;
  ZeroMemory(&ofn, sizeof ofn);
  ofn.lStructSize = sizeof ofn;
  ofn.hwndOwner = owner;
  ofn.lpstrFilter = filter;
  ofn.lpstrFile = buf;
  ofn.nMaxFile = MAX_PATH;
  ofn.Flags = OFN_EXPLORER | OFN_PATHMUSTEXIST | (save ? OFN_OVERWRITEPROMPT : OFN_FILEMUSTEXIST);
  if (save) ofn.lpstrDefExt = "txt";
  return (save ? GetSaveFileNameA(&ofn) : GetOpenFileNameA(&ofn)) ? std::string(buf) : std::string();
}

// fill the accession list from a protein FASTA or a UniProt proteome TSV
static void load_proteins(HWND hwnd, const std::string& path) {
  if (path.empty() || path == g_loaded) return;
  try {
    g_queries = read_queries(path);
  } catch (const std::exception& e) {
    g_queries.clear();
    g_loaded.clear();
    SendMessageA(hAcc, CB_RESETCONTENT, 0, 0);
    SetWindowTextA(hStatus, (std::string("error: ") + e.what()).c_str());
    return;
  }
  g_loaded = path;
  SendMessageA(hAcc, WM_SETREDRAW, FALSE, 0);
  SendMessageA(hAcc, CB_RESETCONTENT, 0, 0);
  SendMessageA(hAcc, CB_INITSTORAGE, (WPARAM)g_queries.size(), 64);
  for (size_t i = 0; i < g_queries.size(); i++) {
    std::string s = g_queries[i].id;
    if (!g_queries[i].desc.empty()) s += "   " + g_queries[i].desc;
    if (s.size() > 110) s.resize(110);
    SendMessageA(hAcc, CB_ADDSTRING, 0, (LPARAM)s.c_str());
  }
  SendMessageA(hAcc, WM_SETREDRAW, TRUE, 0);
  InvalidateRect(hAcc, NULL, TRUE);
  SendMessageA(hAcc, CB_SETCURSEL, 0, 0);
  char st[160];
  snprintf(st, sizeof st, "%d entries - type an accession to jump, or scroll",
           (int)g_queries.size());
  SetWindowTextA(hStatus, st);
  (void)hwnd;
}

static void run_search(HWND hwnd) {
  if (g_busy) return;
  load_proteins(hwnd, get_text(hQ));  // in case the path was typed, not browsed
  std::string dpath = get_text(hD), ev = get_text(hE);
  int sel = (int)SendMessageA(hAcc, CB_GETCURSEL, 0, 0);
  if (g_queries.empty() || sel < 0 || dpath.empty()) {
    MessageBoxA(hwnd, "Choose a protein file, an entry in it, and a genome FASTA.",
                "tblastn-lite", MB_ICONINFORMATION);
    return;
  }
  Opts o;
  o.evalue = ev.empty() ? 1e-5 : atof(ev.c_str());
  o.show_aln = SendMessageA(hAln, BM_GETCHECK, 0, 0) == BST_CHECKED;
  o.show_dna = SendMessageA(hDna, BM_GETCHECK, 0, 0) == BST_CHECKED;
  Rec query = g_queries[(size_t)sel];
  g_busy = true;
  EnableWindow(hRun, FALSE);
  SetWindowTextA(hStatus, ("searching with " + query.id + "...").c_str());
  std::thread([hwnd, dpath, o, query]() {
    std::string* out = new std::string;
    size_t total = 0;
    try {
      Db db = load_db(dpath);
      char buf[512];
      snprintf(buf, sizeof buf, "# tblastn-lite  db: %s (%d seqs, %lld bp)\n", dpath.c_str(),
               (int)db.recs.size(), db.len);
      *out = buf;
      if (!query.desc.empty()) *out += "# " + query.id + "  " + query.desc + "\n";
      std::vector<signed char> q = encode_protein(query.seq);
      std::vector<Hit> hits = search(q, db, o, 0);
      total = hits.size();
      *out += format_hits(query.id, q, db, hits, o);
    } catch (const std::exception& e) {
      *out = std::string("error: ") + e.what() + "\n";
    }
    PostMessageA(hwnd, WM_DONE, (WPARAM)total, (LPARAM)out);
  }).detach();
}

static void layout(HWND hwnd) {
  RECT r;
  GetClientRect(hwnd, &r);
  int W = r.right, pad = 8, lbl = 96, btn = 78, row = 26, y = pad;
  int wide = W - pad * 3 - lbl - btn;
  MoveWindow(hL1, pad, y + 3, lbl - 4, 18, TRUE);
  MoveWindow(hQ, pad + lbl, y, wide, 22, TRUE);
  MoveWindow(hPick, W - pad - btn, y, btn, 22, TRUE);
  y += row;
  MoveWindow(hL2, pad, y + 3, lbl - 4, 18, TRUE);
  MoveWindow(hAcc, pad + lbl, y, W - pad * 2 - lbl, 320, TRUE);  // height = dropdown extent
  y += row;
  MoveWindow(hL3, pad, y + 3, lbl - 4, 18, TRUE);
  MoveWindow(hD, pad + lbl, y, wide, 22, TRUE);
  MoveWindow(hDPick, W - pad - btn, y, btn, 22, TRUE);
  y += row;
  MoveWindow(hL4, pad, y + 3, lbl - 4, 18, TRUE);
  MoveWindow(hE, pad + lbl, y, 80, 22, TRUE);
  MoveWindow(hAln, pad + lbl + 90, y + 2, 100, 20, TRUE);
  MoveWindow(hDna, pad + lbl + 195, y + 2, 110, 20, TRUE);
  MoveWindow(hSave, W - pad * 2 - btn * 2, y, btn, 22, TRUE);
  MoveWindow(hRun, W - pad - btn, y, btn, 22, TRUE);
  y += row + 2;
  MoveWindow(hStatus, pad, r.bottom - 22, W - pad * 2, 18, TRUE);
  MoveWindow(hOut, pad, y, W - pad * 2, r.bottom - y - 26, TRUE);
}

static HWND mk(HWND p, const char* cls, const char* txt, DWORD style, int id) {
  bool sunken = strcmp(cls, "EDIT") == 0;
  HWND h = CreateWindowExA(sunken ? WS_EX_CLIENTEDGE : 0, cls, txt, WS_CHILD | WS_VISIBLE | style,
                           0, 0, 10, 10, p, (HMENU)(INT_PTR)id, NULL, NULL);
  SendMessageA(h, WM_SETFONT, (WPARAM)hUi, TRUE);
  return h;
}

static LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp) {
  switch (msg) {
    case WM_CREATE: {
      hUi = (HFONT)GetStockObject(DEFAULT_GUI_FONT);
      hMono = CreateFontA(-13, 0, 0, 0, FW_NORMAL, 0, 0, 0, ANSI_CHARSET, 0, 0, 0,
                          FIXED_PITCH | FF_MODERN, "Consolas");
      hL1 = mk(hwnd, "STATIC", "Protein file", 0, IDC_L1);
      hQ = mk(hwnd, "EDIT", "", ES_AUTOHSCROLL, IDC_QPATH);
      hPick = mk(hwnd, "BUTTON", "Browse...", 0, IDC_QPICK);
      hL2 = mk(hwnd, "STATIC", "Accession", 0, IDC_L2);
      hAcc = mk(hwnd, "COMBOBOX", "", CBS_DROPDOWNLIST | WS_VSCROLL, IDC_ACC);
      hL3 = mk(hwnd, "STATIC", "Genome FASTA", 0, IDC_L3);
      hD = mk(hwnd, "EDIT", "", ES_AUTOHSCROLL, IDC_DPATH);
      hDPick = mk(hwnd, "BUTTON", "Browse...", 0, IDC_DPICK);
      hL4 = mk(hwnd, "STATIC", "E-value", 0, IDC_L4);
      hE = mk(hwnd, "EDIT", "1e-5", ES_AUTOHSCROLL, IDC_EVAL);
      hAln = mk(hwnd, "BUTTON", "alignments", BS_AUTOCHECKBOX, IDC_ALN);
      hDna = mk(hwnd, "BUTTON", "matching DNA", BS_AUTOCHECKBOX, IDC_DNA);
      hSave = mk(hwnd, "BUTTON", "Save...", 0, IDC_SAVE);
      hRun = mk(hwnd, "BUTTON", "Search", BS_DEFPUSHBUTTON, IDC_RUN);
      hStatus = mk(hwnd, "STATIC", "protein FASTA or UniProt proteome TSV - drag files onto "
                                   "this window, or Browse", 0, IDC_STATUS);
      hOut = mk(hwnd, "EDIT", "",
                ES_MULTILINE | ES_READONLY | ES_AUTOHSCROLL | WS_VSCROLL | WS_HSCROLL, IDC_OUT);
      SendMessageA(hOut, WM_SETFONT, (WPARAM)hMono, TRUE);
      SendMessageA(hAcc, WM_SETFONT, (WPARAM)hMono, TRUE);
      SendMessageA(hOut, EM_SETLIMITTEXT, 0, 0);
      SendMessageA(hAln, BM_SETCHECK, BST_CHECKED, 0);
      DragAcceptFiles(hwnd, TRUE);
      layout(hwnd);
      return 0;
    }
    case WM_SIZE: layout(hwnd); return 0;
    case WM_GETMINMAXINFO: ((MINMAXINFO*)lp)->ptMinTrackSize = POINT{560, 340}; return 0;
    case WM_DROPFILES: {
      HDROP drop = (HDROP)wp;
      char path[MAX_PATH];
      if (DragQueryFileA(drop, 0, path, MAX_PATH)) {
        POINT pt;
        DragQueryPoint(drop, &pt);
        HWND target = ChildWindowFromPoint(hwnd, pt);
        if (target != hQ && target != hD) target = get_text(hQ).empty() ? hQ : hD;
        SetWindowTextA(target, path);
        if (target == hQ) load_proteins(hwnd, path);
      }
      DragFinish(drop);
      return 0;
    }
    case WM_COMMAND:
      if (LOWORD(wp) == IDC_QPATH && HIWORD(wp) == EN_KILLFOCUS) {
        load_proteins(hwnd, get_text(hQ));  // path typed or pasted by hand
        return 0;
      }
      switch (LOWORD(wp)) {
        case IDC_QPICK: {
          std::string p = pick_file(hwnd,
              "Protein FASTA or proteome TSV\0*.fa;*.faa;*.fasta;*.fas;*.pep;*.tsv;*.tab;*.txt\0"
              "All\0*.*\0", false);
          if (p.empty()) return 0;
          SetWindowTextA(hQ, p.c_str());
          load_proteins(hwnd, p);
          return 0;
        }
        case IDC_DPICK: {
          std::string p = pick_file(hwnd,
              "Nucleotide FASTA\0*.fa;*.fna;*.fasta;*.fas;*.ffn;*.txt\0All\0*.*\0", false);
          if (!p.empty()) SetWindowTextA(hD, p.c_str());
          return 0;
        }
        case IDC_RUN: run_search(hwnd); return 0;
        case IDC_SAVE: {
          if (g_text.empty()) return 0;
          std::string p = pick_file(hwnd, "Text\0*.txt\0FASTA\0*.fna\0All\0*.*\0", true);
          if (p.empty()) return 0;
          std::ofstream f(p.c_str(), std::ios::binary);
          f << g_text;
          SetWindowTextA(hStatus, f ? ("saved to " + p).c_str() : "error: could not save");
          return 0;
        }
      }
      return 0;
    case WM_DONE: {
      std::string* res = (std::string*)lp;
      g_text = *res;
      delete res;
      // an EDIT control chokes on megabytes; show the head, Save writes all of it
      const size_t CAP = 1u << 21;
      std::string shown = g_text.size() > CAP
          ? g_text.substr(0, CAP) + "\n[... truncated for display - use Save... for the full result]\n"
          : g_text;
      SetWindowTextA(hOut, to_crlf(shown).c_str());
      char st[128];
      snprintf(st, sizeof st, "done - %d hit%s", (int)wp, wp == 1 ? "" : "s");
      SetWindowTextA(hStatus, st);
      EnableWindow(hRun, TRUE);
      g_busy = false;
      return 0;
    }
    case WM_DESTROY: PostQuitMessage(0); return 0;
  }
  return DefWindowProcA(hwnd, msg, wp, lp);
}

int WINAPI WinMain(HINSTANCE hInst, HINSTANCE, LPSTR, int show) {
  init_tables();
  WNDCLASSA wc;
  ZeroMemory(&wc, sizeof wc);
  wc.lpfnWndProc = WndProc;
  wc.hInstance = hInst;
  wc.hCursor = LoadCursor(NULL, IDC_ARROW);
  wc.hbrBackground = (HBRUSH)(COLOR_BTNFACE + 1);
  wc.lpszClassName = "tblastn_lite";
  wc.hIcon = LoadIcon(NULL, IDI_APPLICATION);
  RegisterClassA(&wc);
  HWND hwnd = CreateWindowExA(WS_EX_ACCEPTFILES, "tblastn_lite", "tblastn-lite",
                              WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT, 940, 640, NULL,
                              NULL, hInst, NULL);
  if (!hwnd) return 1;
  ShowWindow(hwnd, show);
  MSG msg;
  while (GetMessageA(&msg, NULL, 0, 0) > 0) {
    if (!IsDialogMessageA(hwnd, &msg)) {  // tab between fields, Enter = Search
      TranslateMessage(&msg);
      DispatchMessageA(&msg);
    }
  }
  return 0;
}

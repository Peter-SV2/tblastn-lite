// tblastn-lite - plain Win32 GUI.  No toolkit, no runtime: pick two FASTA files, hit Search.
// Build: g++ -O3 -std=c++17 -static -mwindows -o tblastn_gui tblastn_gui.cpp -lcomdlg32 -lshell32
#include <windows.h>
#include <commdlg.h>
#include <shellapi.h>

#include "tblastn_core.h"

enum {
  IDC_QPATH = 101, IDC_QPICK, IDC_DPATH, IDC_DPICK, IDC_EVAL,
  IDC_ALN, IDC_DNA, IDC_RUN, IDC_SAVE, IDC_OUT, IDC_STATUS
};
static const UINT WM_DONE = WM_APP + 1;  // lParam = new std::string* (result text)

static HWND hQ, hD, hE, hAln, hDna, hRun, hSave, hOut, hStatus;
static HFONT hMono, hUi;
static std::string g_text;     // full result, unabridged, for Save
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

static void run_search(HWND hwnd) {
  if (g_busy) return;
  std::string qpath = get_text(hQ), dpath = get_text(hD), ev = get_text(hE);
  if (qpath.empty() || dpath.empty()) {
    MessageBoxA(hwnd, "Choose a protein FASTA and a genome FASTA first.", "tblastn-lite",
                MB_ICONINFORMATION);
    return;
  }
  Opts o;
  o.evalue = ev.empty() ? 1e-5 : atof(ev.c_str());
  o.show_aln = SendMessageA(hAln, BM_GETCHECK, 0, 0) == BST_CHECKED;
  o.show_dna = SendMessageA(hDna, BM_GETCHECK, 0, 0) == BST_CHECKED;
  g_busy = true;
  EnableWindow(hRun, FALSE);
  SetWindowTextA(hStatus, "searching...");
  std::thread([hwnd, qpath, dpath, o]() {
    std::string* out = new std::string;
    size_t total = 0;
    try {
      Db db = load_db(dpath);
      std::vector<Rec> queries = read_fasta(qpath);
      char buf[512];
      snprintf(buf, sizeof buf, "# tblastn-lite  db: %s (%d seqs, %lld bp)\n", dpath.c_str(),
               (int)db.recs.size(), db.len);
      *out = buf;
      for (size_t i = 0; i < queries.size(); i++) {
        std::vector<signed char> q = encode_protein(queries[i].seq);
        if (q.size() < 3) continue;
        std::vector<Hit> hits = search(q, db, o, 0);
        total += hits.size();
        *out += format_hits(queries[i].id, q, db, hits, o);
      }
    } catch (const std::exception& e) {
      *out = std::string("error: ") + e.what() + "\n";
    }
    PostMessageA(hwnd, WM_DONE, (WPARAM)total, (LPARAM)out);
  }).detach();
}

static void layout(HWND hwnd) {
  RECT r;
  GetClientRect(hwnd, &r);
  int W = r.right, pad = 8, lbl = 92, btn = 78, row = 26, y = pad;
  MoveWindow(GetDlgItem(hwnd, IDC_QPATH), pad + lbl, y, W - pad * 3 - lbl - btn, 22, TRUE);
  MoveWindow(GetDlgItem(hwnd, IDC_QPICK), W - pad - btn, y, btn, 22, TRUE);
  y += row;
  MoveWindow(GetDlgItem(hwnd, IDC_DPATH), pad + lbl, y, W - pad * 3 - lbl - btn, 22, TRUE);
  MoveWindow(GetDlgItem(hwnd, IDC_DPICK), W - pad - btn, y, btn, 22, TRUE);
  y += row;
  MoveWindow(hE, pad + lbl, y, 80, 22, TRUE);
  MoveWindow(hAln, pad + lbl + 90, y + 2, 100, 20, TRUE);
  MoveWindow(hDna, pad + lbl + 195, y + 2, 110, 20, TRUE);
  MoveWindow(hRun, W - pad - btn, y, btn, 22, TRUE);
  MoveWindow(hSave, W - pad * 2 - btn * 2, y, btn, 22, TRUE);
  y += row + 2;
  MoveWindow(hStatus, pad, r.bottom - 22, W - pad * 2, 18, TRUE);
  MoveWindow(hOut, pad, y, W - pad * 2, r.bottom - y - 26, TRUE);
}

static HWND mk(HWND p, const char* cls, const char* txt, DWORD style, int id) {
  HWND h = CreateWindowExA(strcmp(cls, "EDIT") == 0 ? WS_EX_CLIENTEDGE : 0, cls, txt,
                           WS_CHILD | WS_VISIBLE | style, 0, 0, 10, 10, p, (HMENU)(INT_PTR)id,
                           NULL, NULL);
  SendMessageA(h, WM_SETFONT, (WPARAM)hUi, TRUE);
  return h;
}

static LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp) {
  switch (msg) {
    case WM_CREATE: {
      hUi = (HFONT)GetStockObject(DEFAULT_GUI_FONT);
      hMono = CreateFontA(-13, 0, 0, 0, FW_NORMAL, 0, 0, 0, ANSI_CHARSET, 0, 0, 0,
                          FIXED_PITCH | FF_MODERN, "Consolas");
      CreateWindowExA(0, "STATIC", "Protein FASTA", WS_CHILD | WS_VISIBLE, 8, 11, 88, 18, hwnd,
                      NULL, NULL, NULL);
      CreateWindowExA(0, "STATIC", "Genome FASTA", WS_CHILD | WS_VISIBLE, 8, 37, 88, 18, hwnd,
                      NULL, NULL, NULL);
      CreateWindowExA(0, "STATIC", "E-value", WS_CHILD | WS_VISIBLE, 8, 63, 88, 18, hwnd, NULL,
                      NULL, NULL);
      for (HWND c = GetWindow(hwnd, GW_CHILD); c; c = GetWindow(c, GW_HWNDNEXT))
        SendMessageA(c, WM_SETFONT, (WPARAM)hUi, TRUE);
      hQ = mk(hwnd, "EDIT", "", ES_AUTOHSCROLL, IDC_QPATH);
      mk(hwnd, "BUTTON", "Browse...", 0, IDC_QPICK);
      hD = mk(hwnd, "EDIT", "", ES_AUTOHSCROLL, IDC_DPATH);
      mk(hwnd, "BUTTON", "Browse...", 0, IDC_DPICK);
      hE = mk(hwnd, "EDIT", "1e-5", ES_AUTOHSCROLL, IDC_EVAL);
      hAln = mk(hwnd, "BUTTON", "alignments", BS_AUTOCHECKBOX, IDC_ALN);
      hDna = mk(hwnd, "BUTTON", "matching DNA", BS_AUTOCHECKBOX, IDC_DNA);
      hSave = mk(hwnd, "BUTTON", "Save...", 0, IDC_SAVE);
      hRun = mk(hwnd, "BUTTON", "Search", BS_DEFPUSHBUTTON, IDC_RUN);
      hStatus = mk(hwnd, "STATIC", "drag FASTA files onto this window, or Browse", 0, IDC_STATUS);
      hOut = mk(hwnd, "EDIT", "", ES_MULTILINE | ES_READONLY | ES_AUTOHSCROLL | WS_VSCROLL |
                                     WS_HSCROLL, IDC_OUT);
      SendMessageA(hOut, WM_SETFONT, (WPARAM)hMono, TRUE);
      SendMessageA(hOut, EM_SETLIMITTEXT, 0, 0);
      SendMessageA(hAln, BM_SETCHECK, BST_CHECKED, 0);
      DragAcceptFiles(hwnd, TRUE);
      layout(hwnd);
      return 0;
    }
    case WM_SIZE: layout(hwnd); return 0;
    case WM_GETMINMAXINFO: ((MINMAXINFO*)lp)->ptMinTrackSize = POINT{520, 320}; return 0;
    case WM_DROPFILES: {
      HDROP drop = (HDROP)wp;
      char path[MAX_PATH];
      if (DragQueryFileA(drop, 0, path, MAX_PATH)) {
        POINT pt;
        DragQueryPoint(drop, &pt);
        HWND target = ChildWindowFromPoint(hwnd, pt);
        if (target != hQ && target != hD) target = get_text(hQ).empty() ? hQ : hD;
        SetWindowTextA(target, path);
      }
      DragFinish(drop);
      return 0;
    }
    case WM_COMMAND:
      switch (LOWORD(wp)) {
        case IDC_QPICK: {
          std::string p = pick_file(hwnd, "FASTA\0*.fa;*.faa;*.fasta;*.fas;*.pep;*.txt\0All\0*.*\0", false);
          if (!p.empty()) SetWindowTextA(hQ, p.c_str());
          return 0;
        }
        case IDC_DPICK: {
          std::string p = pick_file(hwnd, "FASTA\0*.fa;*.fna;*.fasta;*.fas;*.ffn;*.txt\0All\0*.*\0", false);
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
                              WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT, 900, 620, NULL,
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

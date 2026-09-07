// tblastn-lite - Win32 + Direct3D 11 + Dear ImGui, in Echelle's shape.
//
// WinMain, not main: a console flashing behind a window is what the CLI is for.
//
// It does NOT spin. WaitMessage blocks until there is input, so an idle window
// costs nothing; the one exception is while a search is running, when the
// worker has to be able to hand its result back at the next frame.
#include <d3d11.h>
#include <windows.h>
#include <commdlg.h>
#include <shellapi.h>

#include <atomic>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "backends/imgui_impl_dx11.h"
#include "backends/imgui_impl_win32.h"
#include "imgui.h"
#include "imgui_stdlib.h"
#include "tblastn_core.h"

extern IMGUI_IMPL_API LRESULT ImGui_ImplWin32_WndProcHandler(HWND, UINT, WPARAM, LPARAM);

namespace {

ID3D11Device* g_device = nullptr;
ID3D11DeviceContext* g_ctx = nullptr;
IDXGISwapChain* g_swap = nullptr;
ID3D11RenderTargetView* g_rtv = nullptr;
bool g_occluded = false;
UINT g_resize_w = 0, g_resize_h = 0;

// ------------------------------------------------------------------- state ---
// What a finished search hands back. The genome is in here too: the detail pane
// re-translates the frame to draw an alignment, so the sequences have to
// outlive the worker that read them.
struct Result {
  Db db;
  std::vector<Hit> hits;
  std::vector<signed char> q;
  std::string qid, qdesc, err;
  double seconds = 0;
};

struct App {
  std::string protein_path, genome_path;
  std::string filter, evalue = "1e-5";
  std::vector<Rec> entries;      // the loaded FASTA / proteome TSV
  std::vector<int> shown;        // indices of entries matching `filter`
  int sel_entry = -1;            // index into entries, not into shown
  bool aln = true, dna = false;
  std::string status = "open a protein FASTA or a UniProt proteome TSV";

  std::unique_ptr<Result> res;   // the search that is on screen
  int sel_hit = -1;
  std::string detail;            // the selected hit's alignment, as shown
  int detail_key = -1;           // what `detail` was built from

  std::atomic<bool> busy{false};
  std::mutex mtx;
  std::unique_ptr<Result> incoming;  // handed over by the worker thread
};

App g_app;

std::string lower(std::string s) {
  for (size_t i = 0; i < s.size(); i++) s[i] = (char)tolower((unsigned char)s[i]);
  return s;
}

void refilter(App& a) {
  std::string f = lower(a.filter);
  a.shown.clear();
  for (size_t i = 0; i < a.entries.size(); i++) {
    if (f.empty() || lower(a.entries[i].id + " " + a.entries[i].desc).find(f) != std::string::npos)
      a.shown.push_back((int)i);
  }
}

void load_proteins(App& a, const std::string& path) {
  try {
    a.entries = read_queries(path);
  } catch (const std::exception& e) {
    a.entries.clear();
    a.shown.clear();
    a.sel_entry = -1;
    a.status = std::string("error: ") + e.what();
    return;
  }
  a.protein_path = path;
  a.sel_entry = a.entries.size() == 1 ? 0 : -1;
  a.filter.clear();
  refilter(a);
  a.status = std::to_string(a.entries.size()) + " entries - search by accession, gene or name";
}

std::string pick_file(HWND owner, const char* filter, bool save) {
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

// A dropped file goes to the field it belongs in: a TSV is a proteome, and a
// FASTA is sorted by what its first sequence is made of. Asking the user which
// box to aim at is the kind of question the file already answers.
bool looks_nucleotide(const std::string& path) {
  std::ifstream in(path);
  std::string line;
  if (!std::getline(in, line) || line.empty() || line[0] != '>') return false;
  int acgtn = 0, total = 0;
  while (std::getline(in, line) && total < 400) {
    if (!line.empty() && line[0] == '>') break;
    for (size_t i = 0; i < line.size(); i++) {
      if (isspace((unsigned char)line[i])) continue;
      total++;
      if (strchr("ACGTNUacgtnu", line[i])) acgtn++;
    }
  }
  return total > 0 && acgtn * 100 / total >= 85;
}

void start_search(App& a) {
  if (a.busy.load()) return;
  if (a.sel_entry < 0 || a.genome_path.empty()) {
    a.status = "pick an entry and a genome FASTA first";
    return;
  }
  Rec query = a.entries[(size_t)a.sel_entry];
  Opts o;
  o.evalue = a.evalue.empty() ? 1e-5 : atof(a.evalue.c_str());
  std::string gpath = a.genome_path;
  a.busy.store(true);
  a.status = "searching with " + query.id + "...";
  std::thread([&a, query, gpath, o]() {
    std::unique_ptr<Result> r(new Result);
    r->qid = query.id;
    r->qdesc = query.desc;
    LARGE_INTEGER f, t0, t1;
    QueryPerformanceFrequency(&f);
    QueryPerformanceCounter(&t0);
    try {
      r->db = load_db(gpath);
      r->q = encode_protein(query.seq);
      r->hits = search(r->q, r->db, o, 0);
    } catch (const std::exception& e) {
      r->err = e.what();
    }
    QueryPerformanceCounter(&t1);
    r->seconds = double(t1.QuadPart - t0.QuadPart) / double(f.QuadPart);
    std::lock_guard<std::mutex> g(a.mtx);
    a.incoming = std::move(r);
    a.busy.store(false);
  }).detach();
}

// The whole thing on screen as text: what Save writes and what Copy puts on
// the clipboard, so the two cannot disagree.
std::string full_report(const App& a) {
  if (!a.res) return "";
  Opts o;
  o.show_aln = a.aln;
  o.show_dna = a.dna;
  std::string out = "# tblastn-lite  db: " + a.genome_path + "\n";
  if (!a.res->qdesc.empty()) out += "# " + a.res->qid + "  " + a.res->qdesc + "\n";
  return out + format_hits(a.res->qid, a.res->q, a.res->db, a.res->hits, o);
}

// --------------------------------------------------------------------- draw ---
void draw(App& a, HWND hwnd) {
  {
    std::lock_guard<std::mutex> g(a.mtx);
    if (a.incoming) {
      a.res = std::move(a.incoming);
      a.sel_hit = a.res->hits.empty() ? -1 : 0;
      a.detail_key = -1;  // the cached alignment belongs to the old search
      char buf[160];
      if (!a.res->err.empty()) snprintf(buf, sizeof buf, "error: %s", a.res->err.c_str());
      else snprintf(buf, sizeof buf, "%d hits in %.2f s", (int)a.res->hits.size(), a.res->seconds);
      a.status = buf;
    }
  }

  const ImGuiViewport* vp = ImGui::GetMainViewport();
  ImGui::SetNextWindowPos(vp->WorkPos);
  ImGui::SetNextWindowSize(vp->WorkSize);
  ImGui::Begin("tblastn-lite", nullptr,
               ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove |
                   ImGuiWindowFlags_NoBringToFrontOnFocus);

  // --- the rail -------------------------------------------------------------
  ImGui::BeginChild("rail", ImVec2(430, 0), ImGuiChildFlags_Border);

  ImGui::TextUnformatted("Protein  (FASTA or UniProt proteome TSV)");
  ImGui::SetNextItemWidth(-90);
  if (ImGui::InputText("##ppath", &a.protein_path, ImGuiInputTextFlags_EnterReturnsTrue))
    load_proteins(a, a.protein_path);
  ImGui::SameLine();
  if (ImGui::Button("Open...##p", ImVec2(-1, 0))) {
    std::string p = pick_file(hwnd,
        "Protein FASTA or proteome TSV\0*.fa;*.faa;*.fasta;*.fas;*.pep;*.tsv;*.tab;*.txt\0"
        "All\0*.*\0", false);
    if (!p.empty()) load_proteins(a, p);
  }

  ImGui::Spacing();
  ImGui::TextUnformatted("Accession");
  ImGui::SetNextItemWidth(-1);
  if (ImGui::InputTextWithHint("##filter", "search accession, gene or description", &a.filter))
    refilter(a);
  if (!a.entries.empty()) {
    ImGui::TextDisabled("%d of %d entries", (int)a.shown.size(), (int)a.entries.size());
    float h = ImGui::GetTextLineHeightWithSpacing() * 12;
    ImGui::BeginChild("entries", ImVec2(0, h), ImGuiChildFlags_Border);
    ImGuiListClipper clip;  // 4,000-odd rows: only the visible ones are built
    clip.Begin((int)a.shown.size());
    while (clip.Step())
      for (int i = clip.DisplayStart; i < clip.DisplayEnd; i++) {
        int e = a.shown[(size_t)i];
        const Rec& r = a.entries[(size_t)e];
        std::string label = r.id + "   " + r.desc;
        if (ImGui::Selectable(label.c_str(), a.sel_entry == e)) a.sel_entry = e;
        if (ImGui::IsItemHovered() && !r.desc.empty()) ImGui::SetTooltip("%s", r.desc.c_str());
      }
    ImGui::EndChild();
  }
  if (a.sel_entry >= 0) {
    const Rec& r = a.entries[(size_t)a.sel_entry];
    ImGui::Text("query: %s  %d aa", r.id.c_str(), (int)r.seq.size());
  } else {
    ImGui::TextDisabled("no entry selected");
  }

  ImGui::Spacing();
  ImGui::TextUnformatted("Genome  (nucleotide FASTA)");
  ImGui::SetNextItemWidth(-90);
  ImGui::InputText("##gpath", &a.genome_path);
  ImGui::SameLine();
  if (ImGui::Button("Open...##g", ImVec2(-1, 0))) {
    std::string p = pick_file(hwnd,
        "Nucleotide FASTA\0*.fa;*.fna;*.fasta;*.fas;*.ffn;*.txt\0All\0*.*\0", false);
    if (!p.empty()) a.genome_path = p;
  }

  ImGui::Spacing();
  ImGui::SetNextItemWidth(110);
  ImGui::InputText("E-value", &a.evalue);
  ImGui::SameLine();
  ImGui::Checkbox("alignment", &a.aln);
  ImGui::SameLine();
  ImGui::Checkbox("DNA", &a.dna);

  ImGui::Spacing();
  ImGui::BeginDisabled(a.busy.load());
  if (ImGui::Button("Search", ImVec2(120, 0))) start_search(a);
  ImGui::EndDisabled();
  ImGui::SameLine();
  ImGui::BeginDisabled(!a.res || a.res->hits.empty());
  if (ImGui::Button("Copy report", ImVec2(120, 0))) {
    std::string r = full_report(a);
    ImGui::SetClipboardText(r.c_str());
    a.status = "copied " + std::to_string(r.size()) + " bytes to the clipboard";
  }
  ImGui::SameLine();
  if (ImGui::Button("Save report...", ImVec2(-1, 0))) {
    std::string p = pick_file(hwnd, "Text\0*.txt\0FASTA\0*.fna\0All\0*.*\0", true);
    if (!p.empty()) {
      std::ofstream f(p.c_str(), std::ios::binary);
      f << full_report(a);
      a.status = f ? "wrote " + p : "error: could not write " + p;
    }
  }
  ImGui::EndDisabled();

  ImGui::Spacing();
  ImGui::Separator();
  ImGui::TextWrapped("%s", a.status.c_str());
  ImGui::EndChild();

  // --- hits, and the alignment under them -----------------------------------
  ImGui::SameLine();
  ImGui::BeginChild("right", ImVec2(0, 0));
  if (a.res && !a.res->err.empty()) {
    ImGui::TextWrapped("%s", a.res->err.c_str());
  } else if (a.res && a.res->hits.empty()) {
    ImGui::TextDisabled("no hits at this E-value");
  } else if (a.res) {
    const Result& R = *a.res;
    ImGui::Text("%s%s%s   vs   %s", R.qid.c_str(), R.qdesc.empty() ? "" : "  ",
                R.qdesc.c_str(), R.db.recs[0].id.c_str());
    float table_h = ImGui::GetContentRegionAvail().y * (a.aln || a.dna ? 0.42f : 1.0f);
    if (ImGui::BeginTable("hits", 10,
                          ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders |
                              ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingFixedFit,
                          ImVec2(0, table_h))) {
      ImGui::TableSetupScrollFreeze(0, 1);
      const char* cols[] = {"subject", "%id",   "len",    "mism", "frame",
                            "qstart",  "qend",  "sstart", "send", "E-value"};
      for (int c = 0; c < 10; c++) ImGui::TableSetupColumn(cols[c]);
      ImGui::TableHeadersRow();
      for (size_t i = 0; i < R.hits.size(); i++) {
        const Hit& h = R.hits[i];
        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        char lbl[128];
        snprintf(lbl, sizeof lbl, "%s##%d", R.db.recs[h.seqidx].id.c_str(), (int)i);
        if (ImGui::Selectable(lbl, a.sel_hit == (int)i, ImGuiSelectableFlags_SpanAllColumns))
          a.sel_hit = (int)i;
        ImGui::TableNextColumn(); ImGui::Text("%.1f", 100.0 * h.ident / h.len);
        ImGui::TableNextColumn(); ImGui::Text("%d", h.len);
        ImGui::TableNextColumn(); ImGui::Text("%d", h.len - h.ident);
        ImGui::TableNextColumn(); ImGui::Text("%+d", h.frame);
        ImGui::TableNextColumn(); ImGui::Text("%d", h.qs + 1);
        ImGui::TableNextColumn(); ImGui::Text("%d", h.qe + 1);
        ImGui::TableNextColumn(); ImGui::Text("%lld", h.ss);
        ImGui::TableNextColumn(); ImGui::Text("%lld", h.se);
        ImGui::TableNextColumn(); ImGui::Text("%.2g", h.evalue);
      }
      ImGui::EndTable();
    }
    if ((a.aln || a.dna) && a.sel_hit >= 0 && a.sel_hit < (int)R.hits.size()) {
      // Rebuilt from the hit whenever the selection or the switches move, so
      // there is no second copy of an alignment to fall out of step with its
      // row -- but not every frame, because it is also the text being selected
      // in and a rebuild under the cursor would fight the selection.
      int key = a.sel_hit * 4 + (a.aln ? 2 : 0) + (a.dna ? 1 : 0);
      if (key != a.detail_key) {
        a.detail = format_hit_detail(R.q, R.db, R.hits[(size_t)a.sel_hit], a.aln, a.dna);
        a.detail_key = key;
      }
      // A read-only InputTextMultiline, not TextUnformatted: this is output
      // people paste into a notebook, so it has to be selectable with the
      // mouse and answer Ctrl+A / Ctrl+C like any other text box.
      ImGui::InputTextMultiline("##detail", &a.detail, ImVec2(-1, -1),
                                ImGuiInputTextFlags_ReadOnly);
    }
  } else {
    ImGui::TextDisabled("Pick a protein and a genome, then Search.");
  }
  ImGui::EndChild();
  ImGui::End();
}

// ------------------------------------------------------------------- device ---
void create_rtv() {
  ID3D11Texture2D* back = nullptr;
  g_swap->GetBuffer(0, IID_PPV_ARGS(&back));
  if (back) {
    g_device->CreateRenderTargetView(back, nullptr, &g_rtv);
    back->Release();
  }
}

void release_rtv() {
  if (g_rtv) { g_rtv->Release(); g_rtv = nullptr; }
}

bool create_device(HWND hwnd) {
  DXGI_SWAP_CHAIN_DESC sd;
  ZeroMemory(&sd, sizeof sd);
  sd.BufferCount = 2;
  sd.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
  sd.BufferDesc.RefreshRate.Numerator = 60;
  sd.BufferDesc.RefreshRate.Denominator = 1;
  sd.Flags = DXGI_SWAP_CHAIN_FLAG_ALLOW_MODE_SWITCH;
  sd.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
  sd.OutputWindow = hwnd;
  sd.SampleDesc.Count = 1;
  sd.Windowed = TRUE;
  sd.SwapEffect = DXGI_SWAP_EFFECT_DISCARD;
  const D3D_FEATURE_LEVEL want[] = {D3D_FEATURE_LEVEL_11_0, D3D_FEATURE_LEVEL_10_0};
  D3D_FEATURE_LEVEL got;
  HRESULT hr = D3D11CreateDeviceAndSwapChain(nullptr, D3D_DRIVER_TYPE_HARDWARE, nullptr, 0, want,
                                             2, D3D11_SDK_VERSION, &sd, &g_swap, &g_device, &got,
                                             &g_ctx);
  // No GPU, or an RDP session, means no hardware driver. WARP is the software
  // rasteriser and is plenty for this.
  if (hr == DXGI_ERROR_UNSUPPORTED)
    hr = D3D11CreateDeviceAndSwapChain(nullptr, D3D_DRIVER_TYPE_WARP, nullptr, 0, want, 2,
                                       D3D11_SDK_VERSION, &sd, &g_swap, &g_device, &got, &g_ctx);
  if (FAILED(hr)) return false;
  create_rtv();
  return true;
}

void destroy_device() {
  release_rtv();
  if (g_swap) { g_swap->Release(); g_swap = nullptr; }
  if (g_ctx) { g_ctx->Release(); g_ctx = nullptr; }
  if (g_device) { g_device->Release(); g_device = nullptr; }
}

LRESULT WINAPI wndproc(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp) {
  if (ImGui_ImplWin32_WndProcHandler(hwnd, msg, wp, lp)) return true;
  switch (msg) {
    case WM_SIZE:
      if (wp != SIZE_MINIMIZED) { g_resize_w = LOWORD(lp); g_resize_h = HIWORD(lp); }
      return 0;
    case WM_DROPFILES: {
      HDROP drop = (HDROP)wp;
      UINT n = DragQueryFileA(drop, 0xFFFFFFFF, nullptr, 0);
      char path[MAX_PATH];
      for (UINT i = 0; i < n; i++) {
        if (!DragQueryFileA(drop, i, path, MAX_PATH)) continue;
        if (looks_nucleotide(path)) g_app.genome_path = path;
        else load_proteins(g_app, path);
      }
      DragFinish(drop);
      return 0;
    }
    case WM_SYSCOMMAND:
      if ((wp & 0xfff0) == SC_KEYMENU) return 0;  // no Alt menu to open
      return DefWindowProcA(hwnd, msg, wp, lp);
    case WM_DESTROY: PostQuitMessage(0); return 0;
  }
  return DefWindowProcA(hwnd, msg, wp, lp);
}

}  // namespace

int WINAPI WinMain(HINSTANCE inst, HINSTANCE, LPSTR cmdline, int) {
  init_tables();
  WNDCLASSEXA wc;
  ZeroMemory(&wc, sizeof wc);
  wc.cbSize = sizeof wc;
  wc.style = CS_CLASSDC;
  wc.lpfnWndProc = wndproc;
  wc.hInstance = inst;
  wc.hIcon = LoadIcon(nullptr, IDI_APPLICATION);
  wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
  wc.lpszClassName = "tblastn_lite";
  RegisterClassExA(&wc);
  HWND hwnd = CreateWindowA(wc.lpszClassName, "tblastn-lite", WS_OVERLAPPEDWINDOW, 100, 100, 1280,
                            760, nullptr, nullptr, inst, nullptr);
  if (!create_device(hwnd)) {
    destroy_device();
    UnregisterClassA(wc.lpszClassName, inst);
    MessageBoxA(nullptr, "Could not create a Direct3D 11 device.", "tblastn-lite", MB_ICONERROR);
    return 1;
  }
  DragAcceptFiles(hwnd, TRUE);
  ShowWindow(hwnd, SW_SHOWDEFAULT);
  UpdateWindow(hwnd);

  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO& io = ImGui::GetIO();
  io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
  io.IniFilename = nullptr;  // nothing here is worth writing beside the exe
  ImGui::StyleColorsLight();
  ImGuiStyle& st = ImGui::GetStyle();
  st.WindowRounding = 0.0f;
  st.FrameRounding = 2.0f;
  st.ScrollbarRounding = 2.0f;
  st.WindowPadding = ImVec2(8, 8);
  st.ItemSpacing = ImVec2(8, 6);
  ImGui_ImplWin32_Init(hwnd);
  ImGui_ImplDX11_Init(g_device, g_ctx);

  // Files on the command line, so `tblastn_gui proteome.tsv genome.fna` and
  // dropping files on the exe both work. cmdline here excludes the program
  // name, so unlike wWinMain there is no empty-string trap to sidestep.
  if (cmdline && *cmdline) {
    int argc = 0;
    LPWSTR* argv = CommandLineToArgvW(GetCommandLineW(), &argc);
    for (int i = 1; argv && i < argc; i++) {
      char p[MAX_PATH];
      WideCharToMultiByte(CP_ACP, 0, argv[i], -1, p, MAX_PATH, nullptr, nullptr);
      if (looks_nucleotide(p)) g_app.genome_path = p;
      else load_proteins(g_app, p);
    }
    if (argv) LocalFree(argv);
  }

  bool running = true;
  while (running) {
    MSG msg;
    // Block until something happens -- unless a search is in flight, whose
    // result arrives on a thread that posts no message.
    if (!g_app.busy.load() && !PeekMessageA(&msg, nullptr, 0, 0, PM_NOREMOVE)) WaitMessage();
    while (PeekMessageA(&msg, nullptr, 0, 0, PM_REMOVE)) {
      TranslateMessage(&msg);
      DispatchMessageA(&msg);
      if (msg.message == WM_QUIT) running = false;
    }
    if (!running) break;
    if (g_occluded && g_swap->Present(0, DXGI_PRESENT_TEST) == DXGI_STATUS_OCCLUDED) {
      Sleep(10);
      continue;
    }
    g_occluded = false;
    if (g_resize_w != 0 && g_resize_h != 0) {
      release_rtv();
      g_swap->ResizeBuffers(0, g_resize_w, g_resize_h, DXGI_FORMAT_UNKNOWN, 0);
      g_resize_w = g_resize_h = 0;
      create_rtv();
    }
    ImGui_ImplDX11_NewFrame();
    ImGui_ImplWin32_NewFrame();
    ImGui::NewFrame();
    draw(g_app, hwnd);
    ImGui::Render();
    const float bg[4] = {0.941f, 0.941f, 0.941f, 1.0f};
    g_ctx->OMSetRenderTargets(1, &g_rtv, nullptr);
    g_ctx->ClearRenderTargetView(g_rtv, bg);
    ImGui_ImplDX11_RenderDrawData(ImGui::GetDrawData());
    g_occluded = g_swap->Present(1, 0) == DXGI_STATUS_OCCLUDED;
  }

  ImGui_ImplDX11_Shutdown();
  ImGui_ImplWin32_Shutdown();
  ImGui::DestroyContext();
  destroy_device();
  DestroyWindow(hwnd);
  UnregisterClassA(wc.lpszClassName, inst);
  return 0;
}

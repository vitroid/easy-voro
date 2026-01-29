#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <array>
#include <iomanip>

#include "voro++.hh"

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " Lx Ly Lz < points.csv\n";
        return 1;
    }

    // 箱サイズ
    const double Lx = std::atof(argv[1]);
    const double Ly = std::atof(argv[2]);
    const double Lz = std::atof(argv[3]);

    // テキストから点列を読み込み（x y z または x,y,z）
    std::vector<std::array<double, 3>> points;
    points.reserve(10000);

    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty()) continue;

        // 区切り文字としてカンマと空白の両方を許容する
        for (char &ch : line) {
            if (ch == ',') ch = ' ';
        }

        std::stringstream ss(line);
        double x, y, z;
        if (!(ss >> x >> y >> z)) {
            continue;  // パースに失敗した行はスキップ
        }

        points.push_back({x, y, z});
    }

    const int n = static_cast<int>(points.size());
    if (n == 0) {
        std::cerr << "No points read from stdin.\n";
        return 1;
    }

    using namespace voro;

    // Voro++ のコンテナ
    // container(double ax,double bx,double ay,double by,double az,double bz,
    //           int nx,int ny,int nz,bool x_periodic,bool y_periodic,bool z_periodic,int init_mem);
    const int nx = 6, ny = 6, nz = 6;  // 適当な分割数（必要に応じて調整）
    container con(
        0.0, Lx,
        0.0, Ly,
        0.0, Lz,
        nx, ny, nz,
        true, true, true,  // 三方向とも周期
        8                  // 初期メモリ
    );

    // 粒子を投入 (ID は 0..n-1)
    for (int i = 0; i < n; ++i) {
        const double x = points[i][0];
        const double y = points[i][1];
        const double z = points[i][2];
        con.put(i, x, y, z);
    }

    // Voronoi 解析
    using voro::voronoicell_neighbor;

    c_loop_all cl(con);
    voronoicell_neighbor cell;

    double total_volume = 0.0;
    std::cout << std::setprecision(17);

    if (cl.start()) do {
        if (con.compute_cell(cell, cl)) {
            const int id = cl.pid();  // 粒子 ID
            double x, y, z;
            cl.pos(x, y, z);          // 粒子位置（セルの原点）

            const double vol = cell.volume();
            total_volume += vol;

            // 頂点座標（絶対座標系）
            std::vector<double> vv;
            cell.vertices(x, y, z, vv);  // [x0,y0,z0,x1,y1,z1,...]

            // 面ごとの頂点インデックスと隣接セルID
            std::vector<int> face_orders;
            std::vector<int> face_verts;
            std::vector<int> neigh;
            cell.face_orders(face_orders);
            cell.face_vertices(face_verts);
            cell.neighbors(neigh);

            // JSON 1行出力
            std::cout << "{";

            // id
            std::cout << "\"id\":" << id << ",";

            // original position
            std::cout << "\"position\":[" << x << "," << y << "," << z << "],";

            // volume
            std::cout << "\"volume\":" << vol << ",";

            // vertices: [[x,y,z],...]
            std::cout << "\"vertices\":[";
            const int nv = static_cast<int>(vv.size() / 3);
            for (int vi = 0; vi < nv; ++vi) {
                if (vi > 0) std::cout << ",";
                std::cout << "[" << vv[3 * vi] << "," << vv[3 * vi + 1] << "," << vv[3 * vi + 2] << "]";
            }
            std::cout << "],";

            // faces: [{adjacent_cell: id, vertices: [..]}, ...]
            std::cout << "\"faces\":[";
            std::size_t idx = 0;
            const int nf = static_cast<int>(face_orders.size());
            for (int fi = 0; fi < nf; ++fi) {
                if (fi > 0) std::cout << ",";
                const int order = face_orders[fi];
                std::cout << "{";
                // neighbor id
                int nid = (fi < static_cast<int>(neigh.size())) ? neigh[fi] : -1;
                std::cout << "\"adjacent_cell\":" << nid << ",";
                // vertex indices for this face
                std::cout << "\"vertices\":[";
                for (int k = 0; k < order; ++k) {
                    if (k > 0) std::cout << ",";
                    std::cout << face_verts[idx + k];
                }
                std::cout << "]}";
                idx += order;
            }
            std::cout << "]}\n";
        }
    } while (cl.inc());

    const double box_volume = Lx * Ly * Lz;
    std::cerr << "N = " << n << "\n";
    std::cerr << "box_volume    = " << box_volume << "\n";
    std::cerr << "total_volume  = " << total_volume << "\n";
    std::cerr << "ratio (total/box) = " << (total_volume / box_volume) << "\n";

    return 0;
}
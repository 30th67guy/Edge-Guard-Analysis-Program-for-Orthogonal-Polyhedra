#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/IO/Color.h>
#include <GLFW/glfw3.h>
#include <cstdlib>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <cmath>
#include <random>

#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <chrono>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
using Surface_mesh = CGAL::Surface_mesh<Kernel::Point_3>;
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;

namespace PMP = CGAL::Polygon_mesh_processing;
using Color = CGAL::IO::Color;

// mouse rotation
float rotation_x = 0.0f;
float rotation_y = 0.0f;
double last_cursor_x, last_cursor_y;
bool mouse_pressed = false;

Nef_polyhedron make_cube(Point center, double length, double width, double depth) {
    length /= 2.0;
    width /= 2.0;
    depth /= 2.0;

    std::vector<Point> points = {
        Point(center.x()-length, center.y()-width, center.z()-depth),
        Point(center.x()+length, center.y()-width, center.z()-depth),
        Point(center.x()+length, center.y()+width, center.z()-depth),
        Point(center.x()-length, center.y()+width, center.z()-depth),
        Point(center.x()-length, center.y()-width, center.z()+depth),
        Point(center.x()+length, center.y()-width, center.z()+depth),
        Point(center.x()+length, center.y()+width, center.z()+depth),
        Point(center.x()-length, center.y()+width, center.z()+depth),
    };
    Surface_mesh mesh;
    auto v = [&](int i) { return mesh.add_vertex(points[i]); };
    std::vector<Surface_mesh::Vertex_index> vi(8);
    for (int i = 0; i < 8; ++i) vi[i] = v(i);

    auto add_face = [&](int a, int b, int c, int d) {
        mesh.add_face(vi[a], vi[b], vi[c]);
        mesh.add_face(vi[a], vi[c], vi[d]);
    };

    add_face(0, 1, 2, 3);
    add_face(7, 6, 5, 4);
    add_face(0, 4, 5, 1);
    add_face(1, 5, 6, 2);
    add_face(2, 6, 7, 3);
    add_face(3, 7, 4, 0);

    return Nef_polyhedron(mesh);
}

Nef_polyhedron make_regular_ngon_prism(Point prism_center, int num_base_vertices, double r, double height) {
    if (num_base_vertices < 3) {
        std::cerr << "Error: Prism must have at least 3 base vertices. Requested: " << num_base_vertices << std::endl;
        return Nef_polyhedron();
        // Return an empty Nef_polyhedron for invalid input

    }

    Surface_mesh mesh;
    std::vector<Surface_mesh::Vertex_index> base_vertices_vi(num_base_vertices);
    std::vector<Surface_mesh::Vertex_index> top_vertices_vi(num_base_vertices);

    height /= 2.0;

    // Generate and add all vertices to the mesh
    for (int i = 0; i < num_base_vertices; ++i) {
        // Calculate 2D coordinates for a regular n-gon vertex
        double angle = 2.0 * CGAL_PI * static_cast<double>(i) / static_cast<double>(num_base_vertices);
        double x_offset = r * std::cos(angle);
        double z_offset = r * std::sin(angle);

        // Create 3D points for base and top, then add to mesh
        Point base_pt(prism_center.x() + x_offset, prism_center.y() - height, prism_center.z() + z_offset);
        base_vertices_vi[i] = mesh.add_vertex(base_pt);

        Point top_pt(prism_center.x() + x_offset, prism_center.y() + height, prism_center.z() + z_offset);
        top_vertices_vi[i] = mesh.add_vertex(top_pt);
    }

    // triangulated base faces
    // Vertices are generated counter clockwise ( when viewed from +Z )
    for (int i = 1; i < num_base_vertices - 1; ++i) {
        mesh.add_face(base_vertices_vi[0], base_vertices_vi[i + 1], base_vertices_vi[i]);
    }

    // triangulated top faces
    // For the top face (normal pointing in +Z direction)
    for (int i = 1; i < num_base_vertices - 1; ++i) {
        mesh.add_face(top_vertices_vi[0], top_vertices_vi[i], top_vertices_vi[i + 1]);
    }

    // side faces
    // each quadrilateral side is split into two triangles
    for (int i = 0; i < num_base_vertices; ++i) {
        int next_i = (i + 1) % num_base_vertices; // Next vertex index, wraps around

        // Triangle 1: (base_vi[i], base_vi[next_i], top_vi[i])
        // Triangle 2: (base_vi[next_i], top_vi[next_i], top_vi[i])
        mesh.add_face(base_vertices_vi[i], base_vertices_vi[next_i], top_vertices_vi[i]);
        mesh.add_face(base_vertices_vi[next_i], top_vertices_vi[next_i], top_vertices_vi[i]);
    }

    return Nef_polyhedron(mesh);
}

Surface_mesh rotate_mesh(const Surface_mesh& input_mesh, const Kernel::Vector_3& axis, double angle_degrees) {
    Surface_mesh output_mesh = input_mesh;
    Kernel::Vector_3 n = axis / std::sqrt(CGAL::to_double(axis.squared_length()));
    double angle_radians = angle_degrees * CGAL_PI / 180.0;
    double c = std::cos(angle_radians);
    double s = std::sin(angle_radians);
    double t = 1.0 - c;
    double x = CGAL::to_double(n.x()), y = CGAL::to_double(n.y()), z = CGAL::to_double(n.z());

    // Rodrigues' rotation formula
    CGAL::Aff_transformation_3<Kernel> rotation(
        x * x * t + c,
        x * y * t - z * s,
        x * z * t + y * s,
        y * x * t + z * s,
        y * y * t + c,
        y * z * t - x * s,
        z * x * t - y * s,
        z * y * t + x * s,
        z * z * t + c
    );
    for (auto v : output_mesh.vertices()) {
        output_mesh.point(v) = rotation.transform(output_mesh.point(v));
    }
    return output_mesh;
}

Nef_polyhedron rotate_nef_polyhedron(const Nef_polyhedron& nef, const Kernel::Vector_3& axis, double angle_degrees) {
    Surface_mesh mesh;
    CGAL::convert_nef_polyhedron_to_polygon_mesh(nef, mesh);
    Surface_mesh rotated_mesh = rotate_mesh(mesh, axis, angle_degrees);
    return Nef_polyhedron(rotated_mesh);
}

void set_perspective(float fovy, float aspect, float zNear, float zFar) {
    float f = 1.0f / tanf(fovy * 0.5f * 3.14159265f / 180.0f);
    float m[16] = { 0 };
    m[0] = f / aspect;
    m[5] = f;
    m[10] = (zFar + zNear) / (zNear - zFar);
    m[11] = -1.0f;
    m[14] = (2.0f * zFar * zNear) / (zNear - zFar);
    glMultMatrixf(m);
}

void tag_original_edges(Surface_mesh& mesh) {
    auto edge_map_pair = mesh.add_property_map<Surface_mesh::Edge_index, bool>("e:original", false);
    auto& edge_map = edge_map_pair.first;
    for (auto edge : mesh.edges())
        edge_map[edge] = true;
}

void assign_random_face_colors(Surface_mesh& mesh) {
    using Face_index = Surface_mesh::Face_index;
    auto color_map_pair = mesh.add_property_map<Face_index, Color>("f:color", Color(0, 0, 0));
    auto& color_map = color_map_pair.first;
    for (auto face : mesh.faces())
        color_map[face] = Color(
            static_cast<unsigned char>(rand() % 256),
            static_cast<unsigned char>(rand() % 256),
            static_cast<unsigned char>(rand() % 256)
        );
}

void triangulate_faces_with_colors(Surface_mesh& mesh) {
    tag_original_edges(mesh);
    assign_random_face_colors(mesh);
    using Face_index = Surface_mesh::Face_index;
    using Vertex_index = Surface_mesh::Vertex_index;
    auto color_map_opt = mesh.property_map<Face_index, Color>("f:color");
    if (!color_map_opt.has_value())
        throw std::runtime_error("Missing 'f:color' map.");
    auto& color_map = *color_map_opt;
    std::map<Face_index, Color> face_colors;
    std::vector<std::tuple<Point, Point, Point, Color>> triangles;
    for (auto face : mesh.faces()) {
        Color col = color_map[face];
        std::vector<Point> points;
        for (auto v : CGAL::vertices_around_face(mesh.halfedge(face), mesh))
            points.push_back(mesh.point(v));
        if (points.size() >= 3)
            for (std::size_t i = 1; i + 1 < points.size(); ++i)
                triangles.emplace_back(points[0], points[i], points[i + 1], col);
    }
    std::set<std::pair<Point, Point>> original_edges;
    for (auto edge : mesh.edges()) {
        auto h = mesh.halfedge(edge);
        Point p1 = mesh.point(mesh.source(h));
        Point p2 = mesh.point(mesh.target(h));
        if (p1 > p2) std::swap(p1, p2); // Ensure consistent ordering
        original_edges.insert({p1, p2});
    }
    mesh.clear();
    auto new_color_map_pair = mesh.add_property_map<Face_index, Color>("f:color", Color(0, 0, 0));
    auto& new_color_map = new_color_map_pair.first;
    std::map<Point, Vertex_index> point_map;
    auto add_vertex = [&](const Point& p) -> Vertex_index {
        auto it = point_map.find(p);
        if (it != point_map.end()) return it->second;
        Vertex_index v = mesh.add_vertex(p);
        point_map[p] = v;
        return v;
    };
    for (const auto& [p0, p1, p2, col] : triangles) {
        Vertex_index v0 = add_vertex(p0);
        Vertex_index v1 = add_vertex(p1);
        Vertex_index v2 = add_vertex(p2);
        Face_index f = mesh.add_face(v0, v1, v2);
        if (f != Surface_mesh::null_face()) new_color_map[f] = col;
    }
    auto new_edge_map_pair = mesh.add_property_map<Surface_mesh::Edge_index, bool>("e:original", false);
    auto& new_edge_map = new_edge_map_pair.first;
    for (auto edge : mesh.edges()) {
        auto h = mesh.halfedge(edge);
        Point p1 = mesh.point(mesh.source(h));
        Point p2 = mesh.point(mesh.target(h));
        if (p1 > p2) std::swap(p1, p2);
        if (original_edges.count({p1, p2}) > 0)
            new_edge_map[edge] = true;
    }
}

void draw_mesh(const Surface_mesh& mesh) {
    auto color_map_opt = mesh.property_map<Surface_mesh::Face_index, Color>("f:color");
    if (!color_map_opt.has_value())
        throw std::runtime_error("Face color property map 'f:color' not found.");
    const auto& color_map = *color_map_opt;
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0f, 1.0f);
    glBegin(GL_TRIANGLES);
    for (auto face : mesh.faces()) {
        const Color& color = color_map[face];
        glColor3f(
            color.red()   / 255.0f,
            color.green() / 255.0f,
            color.blue()  / 255.0f
        );
        std::vector<Point> verts;
        for (auto v : CGAL::vertices_around_face(mesh.halfedge(face), mesh))
            verts.push_back(mesh.point(v));
        if (verts.size() == 3)
            for (const Point& p : verts)
                glVertex3f(
                    static_cast<float>(CGAL::to_double(p.x())),
                    static_cast<float>(CGAL::to_double(p.y())),
                    static_cast<float>(CGAL::to_double(p.z()))
                );
    }
    glEnd();
    glDisable(GL_POLYGON_OFFSET_FILL);
    
    auto original_map_opt = mesh.property_map<Surface_mesh::Edge_index, bool>("e:original");
    if (!original_map_opt) return;
    const auto& original_map = *original_map_opt;
    
    auto highlight_map_opt = mesh.property_map<Surface_mesh::Edge_index, bool>("e:highlight");

    // --- パス1: ハイライトされていない辺を黄色・通常の太さで描画 ---
    glLineWidth(0.4f); // 通常の線の太さ
    glColor3f(1.0f, 1.0f, 0.0f); // 黄色
    glBegin(GL_LINES);
    for (auto edge : mesh.edges()) {
        if (!original_map[edge]) continue; // 本物の辺でなければスキップ

        bool is_highlighted = highlight_map_opt && (*highlight_map_opt)[edge];
        if (is_highlighted) continue; // ハイライトされている辺は、このパスではスキップ

        auto h = mesh.halfedge(edge);
        const Point& p1 = mesh.point(mesh.source(h));
        const Point& p2 = mesh.point(mesh.target(h));
        glVertex3f(static_cast<float>(CGAL::to_double(p1.x())), static_cast<float>(CGAL::to_double(p1.y())), static_cast<float>(CGAL::to_double(p1.z())));
        glVertex3f(static_cast<float>(CGAL::to_double(p2.x())), static_cast<float>(CGAL::to_double(p2.y())), static_cast<float>(CGAL::to_double(p2.z())));
    }
    glEnd();

    // --- パス2: ハイライトされている辺を赤色・太い線で描画 ---
    if (highlight_map_opt) { // ハイライトマップが存在する場合のみ
        glLineWidth(5.0f); // 太い線の太さ
        glColor3f(1.0f, 0.0f, 0.0f); // 赤色
        glBegin(GL_LINES);
        for (auto edge : mesh.edges()) {
            if (!original_map[edge]) continue; // 本物の辺でなければスキップ

            bool is_highlighted = (*highlight_map_opt)[edge];
            if (!is_highlighted) continue; // ハイライトされていない辺は、このパスではスキップ

            auto h = mesh.halfedge(edge);
            const Point& p1 = mesh.point(mesh.source(h));
            const Point& p2 = mesh.point(mesh.target(h));
            glVertex3f(static_cast<float>(CGAL::to_double(p1.x())), static_cast<float>(CGAL::to_double(p1.y())), static_cast<float>(CGAL::to_double(p1.z())));
            glVertex3f(static_cast<float>(CGAL::to_double(p2.x())), static_cast<float>(CGAL::to_double(p2.y())), static_cast<float>(CGAL::to_double(p2.z())));
        }
        glEnd();
    }
}

// mouse callback
void cursor_position_callback(GLFWwindow* window, double xpos, double ypos) {
    if (mouse_pressed) {
        float dx = static_cast<float>(xpos - last_cursor_x);
        float dy = static_cast<float>(ypos - last_cursor_y);
        rotation_y += dx * 0.2f;
        rotation_x += dy * 0.2f;
    }
    last_cursor_x = xpos;
    last_cursor_y = ypos;
}

// mouse button callback
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS)
            mouse_pressed = true;
        else if (action == GLFW_RELEASE)
            mouse_pressed = false;
    }
}

std::random_device seed_gen;
std::mt19937 mt_engine(seed_gen());
double get_rand(int min_val, int max_val, int ratio, std::mt19937& engine) {
    if (min_val > max_val) {
        std::swap(min_val, max_val);
    }

    std::uniform_int_distribution<int> distrib(min_val, max_val);
    return static_cast<double>(distrib(engine)) / static_cast<double>(ratio);
}

double xyz(std::mt19937& engine) {
    std::uniform_int_distribution<> distrib(0, 1); // 0か1を生成
    return (distrib(engine) == 0) ? 1.0 : -1.0;
}

struct EdgeAnalysisResult {
    int real_edges = 0;   // Amount of Real Edges
    int reflex_edges = 0; // Amount of Reflex Edges in Real Edges
    int convex_edges = 0; // Amount of Convex Edges in Real Edges

    size_t x_aligned_edges = 0; // X軸に平行な辺の数
    size_t y_aligned_edges = 0; // Y軸に平行な辺の数
    size_t z_aligned_edges = 0; // Z軸に平行な辺の数
};

// Define edge types (Fake, Convex, Reflex)
// 辺のタイプを定義 (偽の辺、凸な辺、凹な辺)
enum class EdgeType { FAKE, CONVEX, REFLEX };

// Define directions of the edges
// 辺の軸方向を定義
enum class EdgeAxis { X, Y, Z, OTHER };

// Define 24 patterns which any vertices would be classified
// P means positive, N means negative
// 頂点の「角」が取りうる24種のパターンを定義
enum class CornerPattern {
    UNDEFINED,
    CONVEX_PX_PY, CONVEX_PX_NY, CONVEX_NX_PY, CONVEX_NX_NY,
    CONVEX_PY_PZ, CONVEX_PY_NZ, CONVEX_NY_PZ, CONVEX_NY_NZ,
    CONVEX_PZ_PX, CONVEX_PZ_NX, CONVEX_NZ_PX, CONVEX_NZ_NX,
    REFLEX_PX_PY, REFLEX_PX_NY, REFLEX_NX_PY, REFLEX_NX_NY,
    REFLEX_PY_PZ, REFLEX_PY_NZ, REFLEX_NY_PZ, REFLEX_NY_NZ,
    REFLEX_PZ_PX, REFLEX_PZ_NX, REFLEX_NZ_PX, REFLEX_NZ_NX
};

// 12種類のターゲットパターンセットを定義・初期化する関数
std::vector<std::set<CornerPattern>> define_target_pattern_sets() {
    std::vector<std::set<CornerPattern>> targets(12);

    // --- XY平面上の4セット ---
    // Set 0: 基準 (0度回転)
    targets[0] = {CornerPattern::CONVEX_PX_NY, CornerPattern::REFLEX_PX_PY, CornerPattern::REFLEX_NX_NY};
    // Set 1: 90度回転 (+X -> +Y, +Y -> -X)
    targets[1] = {CornerPattern::CONVEX_PX_PY, CornerPattern::REFLEX_NX_PY, CornerPattern::REFLEX_PX_NY};
    // Set 2: 180度回転
    targets[2] = {CornerPattern::CONVEX_NX_PY, CornerPattern::REFLEX_NX_NY, CornerPattern::REFLEX_PX_PY};
    // Set 3: 270度回転
    targets[3] = {CornerPattern::CONVEX_NX_NY, CornerPattern::REFLEX_PX_NY, CornerPattern::REFLEX_NX_PY};

    // --- YZ平面上の4セット (X->Y, Y->Z, Z->X の置換) ---
    // Set 4: 基準
    targets[4] = {CornerPattern::CONVEX_PY_NZ, CornerPattern::REFLEX_PY_PZ, CornerPattern::REFLEX_NY_NZ};
    // Set 5: 90度回転
    targets[5] = {CornerPattern::CONVEX_PY_PZ, CornerPattern::REFLEX_NY_PZ, CornerPattern::REFLEX_PY_NZ};
    // Set 6: 180度回転
    targets[6] = {CornerPattern::CONVEX_NY_PZ, CornerPattern::REFLEX_NY_NZ, CornerPattern::REFLEX_PY_PZ};
    // Set 7: 270度回転
    targets[7] = {CornerPattern::CONVEX_NY_NZ, CornerPattern::REFLEX_PY_NZ, CornerPattern::REFLEX_NY_PZ};

    // --- ZX平面上の4セット (X->Z, Y->X, Z->Y の置換) ---
    // Set 8: 基準
    targets[8] = {CornerPattern::CONVEX_PZ_NX, CornerPattern::REFLEX_PZ_PX, CornerPattern::REFLEX_NZ_NX};
    // Set 9: 90度回転
    targets[9] = {CornerPattern::CONVEX_PZ_PX, CornerPattern::REFLEX_NZ_PX, CornerPattern::REFLEX_PZ_NX};
    // Set 10: 180度回転
    targets[10] = {CornerPattern::CONVEX_NZ_PX, CornerPattern::REFLEX_NZ_NX, CornerPattern::REFLEX_PZ_PX};
    // Set 11: 270度回転
    targets[11] = {CornerPattern::CONVEX_NZ_NX, CornerPattern::REFLEX_PZ_NX, CornerPattern::REFLEX_NZ_PX};
    
    return targets;
}

// 最終的な集計を行う関数
std::vector<int> count_pattern_sets(
    const Surface_mesh& mesh,
    const std::map<Surface_mesh::Vertex_index, std::vector<CornerPattern>>& vertex_patterns_map)
{
    // 1. ターゲットとなる12セットを定義
    const auto target_sets = define_target_pattern_sets();

    // 2. 12セットそれぞれのカウンターを準備
    std::vector<int> counts(12, 0);

    // 3. 全ての頂点をループしてパターンをチェック
    for (const auto& pair : vertex_patterns_map) {
        std::cout << "Loop Test 1" << std::endl;
        const auto& vertex_patterns = pair.second;

        // 頂点が3つの角パターンを持っている場合のみ処理
        if (vertex_patterns.size() != 3) {
            std::cout << "Loop Test 2" << std::endl;
            continue;
        }

        // 4. 頂点のパターンを順序不問の `std::set` に変換
        std::set<CornerPattern> vertex_set(vertex_patterns.begin(), vertex_patterns.end());

        // 5. 12個のターゲットセットのいずれかに一致するかをチェック
        for (int i = 0; i < target_sets.size(); ++i) {
            std::cout << "Loop Test 3" << std::endl;
            if (vertex_set == target_sets[i]) {
                // 一致したら、対応するカウンターを増やして、この頂点のチェックは終了
                std::cout << "Loop Test 4" << std::endl;
                counts[i]++;
                break;
            }
        }
    }

    // 6. カウント結果を返す
    return counts;
}

// 法線ベクトルが6つのどの軸方向を向いているかを判定するヘルパー関数
enum class NormalDirection { PX, NX, PY, NY, PZ, NZ, OTHER };
NormalDirection get_normal_direction(const Kernel::Vector_3& normal) {
    // 各成分の絶対値を取得
    double x_abs = std::abs(CGAL::to_double(normal.x()));
    double y_abs = std::abs(CGAL::to_double(normal.y()));
    double z_abs = std::abs(CGAL::to_double(normal.z()));

    // 最も大きい成分がどの軸かを判断する
    if (x_abs > y_abs && x_abs > z_abs) { // X軸方向が支配的
        return (CGAL::to_double(normal.x()) > 0) ? NormalDirection::PX : NormalDirection::NX;
    } else if (y_abs > x_abs && y_abs > z_abs) { // Y軸方向が支配的
        return (CGAL::to_double(normal.y()) > 0) ? NormalDirection::PY : NormalDirection::NY;
    } else if (z_abs > x_abs && z_abs > y_abs) { // Z軸方向が支配的
        return (CGAL::to_double(normal.z()) > 0) ? NormalDirection::PZ : NormalDirection::NZ;
    }
    
    return NormalDirection::OTHER; // どの軸にも平行でない場合（今回は起こらないはず）
}

// 各「本物の辺」を24パターンのいずれかに直接分類する最終版の関数
std::map<Surface_mesh::Edge_index, CornerPattern> classify_edges_into_patterns(const Surface_mesh& mesh) {
    std::map<Surface_mesh::Edge_index, CornerPattern> edge_pattern_map;

    for (auto edge_desc : mesh.edges()) {
        auto hf1 = mesh.halfedge(edge_desc, 0);
        if (mesh.is_border(hf1)) { continue; }

        auto f1 = mesh.face(hf1);
        auto f2 = mesh.face(mesh.opposite(hf1));
        Kernel::Vector_3 n1 = CGAL::Polygon_mesh_processing::compute_face_normal(f1, mesh);
        Kernel::Vector_3 n2 = CGAL::Polygon_mesh_processing::compute_face_normal(f2, mesh);

        // --- 1. 「本物の辺」かどうかを判定 ---
        if (CGAL::cross_product(n1, n2).squared_length() == 0) {
            continue; // 偽の辺はスキップ
        }

        // --- ここから下は「本物の辺」に対する処理 ---

        // --- 2. 辺が Convex か Reflex かを判定 ---
        auto v_source = mesh.source(hf1);
        auto v_other_on_f2 = mesh.target(mesh.next(mesh.opposite(hf1)));
        bool is_reflex = ((mesh.point(v_other_on_f2) - mesh.point(v_source)) * n1 > 0);

        // --- 3. 2つの面の法線方向を取得（新しい頑健版の関数を使用）---
        NormalDirection dir1 = get_normal_direction(n1);
        NormalDirection dir2 = get_normal_direction(n2);
        
        // --- 4. 判定結果を組み合わせて、24パターンのどれかを決定 ---
        CornerPattern found_pattern = CornerPattern::UNDEFINED;
        if ((dir1 == NormalDirection::PX && dir2 == NormalDirection::PY) || (dir1 == NormalDirection::PY && dir2 == NormalDirection::PX)) found_pattern = is_reflex ? CornerPattern::REFLEX_PX_PY : CornerPattern::CONVEX_PX_PY;
        else if ((dir1 == NormalDirection::PX && dir2 == NormalDirection::NY) || (dir1 == NormalDirection::NY && dir2 == NormalDirection::PX)) found_pattern = is_reflex ? CornerPattern::REFLEX_PX_NY : CornerPattern::CONVEX_PX_NY;
        else if ((dir1 == NormalDirection::NX && dir2 == NormalDirection::PY) || (dir1 == NormalDirection::PY && dir2 == NormalDirection::NX)) found_pattern = is_reflex ? CornerPattern::REFLEX_NX_PY : CornerPattern::CONVEX_NX_PY;
        else if ((dir1 == NormalDirection::NX && dir2 == NormalDirection::NY) || (dir1 == NormalDirection::NY && dir2 == NormalDirection::NX)) found_pattern = is_reflex ? CornerPattern::REFLEX_NX_NY : CornerPattern::CONVEX_NX_NY;
        else if ((dir1 == NormalDirection::PY && dir2 == NormalDirection::PZ) || (dir1 == NormalDirection::PZ && dir2 == NormalDirection::PY)) found_pattern = is_reflex ? CornerPattern::REFLEX_PY_PZ : CornerPattern::CONVEX_PY_PZ;
        else if ((dir1 == NormalDirection::PY && dir2 == NormalDirection::NZ) || (dir1 == NormalDirection::NZ && dir2 == NormalDirection::PY)) found_pattern = is_reflex ? CornerPattern::REFLEX_PY_NZ : CornerPattern::CONVEX_PY_NZ;
        else if ((dir1 == NormalDirection::NY && dir2 == NormalDirection::PZ) || (dir1 == NormalDirection::PZ && dir2 == NormalDirection::NY)) found_pattern = is_reflex ? CornerPattern::REFLEX_NY_PZ : CornerPattern::CONVEX_NY_PZ;
        else if ((dir1 == NormalDirection::NY && dir2 == NormalDirection::NZ) || (dir1 == NormalDirection::NZ && dir2 == NormalDirection::NY)) found_pattern = is_reflex ? CornerPattern::REFLEX_NY_NZ : CornerPattern::CONVEX_NY_NZ;
        else if ((dir1 == NormalDirection::PZ && dir2 == NormalDirection::PX) || (dir1 == NormalDirection::PX && dir2 == NormalDirection::PZ)) found_pattern = is_reflex ? CornerPattern::REFLEX_PZ_PX : CornerPattern::CONVEX_PZ_PX;
        else if ((dir1 == NormalDirection::PZ && dir2 == NormalDirection::NX) || (dir1 == NormalDirection::NX && dir2 == NormalDirection::PZ)) found_pattern = is_reflex ? CornerPattern::REFLEX_PZ_NX : CornerPattern::CONVEX_PZ_NX;
        else if ((dir1 == NormalDirection::NZ && dir2 == NormalDirection::PX) || (dir1 == NormalDirection::PX && dir2 == NormalDirection::NZ)) found_pattern = is_reflex ? CornerPattern::REFLEX_NZ_PX : CornerPattern::CONVEX_NZ_PX;
        else if ((dir1 == NormalDirection::NZ && dir2 == NormalDirection::NX) || (dir1 == NormalDirection::NX && dir2 == NormalDirection::NZ)) found_pattern = is_reflex ? CornerPattern::REFLEX_NZ_NX : CornerPattern::CONVEX_NZ_NX;

        edge_pattern_map[edge_desc] = found_pattern;
    }
    return edge_pattern_map;
}

void tag_highlight_edges(
    Surface_mesh& mesh,
    const std::map<Surface_mesh::Edge_index, CornerPattern>& edge_pattern_map,
    const std::set<CornerPattern>& highlight_patterns) 
{
    // ハイライトするかどうかを記録する新しいプロパティマップを作成
    auto highlight_map = mesh.add_property_map<Surface_mesh::Edge_index, bool>("e:highlight", false).first;

    // 分類済みの全ての辺をループ
    for (const auto& pair : edge_pattern_map) {
        auto edge_desc = pair.first;
        auto pattern = pair.second;

        // この辺のパターンが、ハイライト対象のセットに含まれているかチェック
        if (highlight_patterns.count(pattern)) {
            // 含まれていれば、ハイライトマップの値を true にする
            highlight_map[edge_desc] = true;
        }
    }
}

struct settings{
        // -- mode
        // when 1, we check the intersections of generated cuboids with intersection function. (Set as default to 1)
        int mode = 1;

        // -- shape
        // when shape = 0, we generate each normal cuboids which has each coordinates and sizes randomly. (Set as default to 0)
        int shape = 0;

        // -- skip_log
        // when false, it won't show logs for debug. (set as default to false)
        bool skip_log = false;

        // -- amount
        // the amount how many cuboids would be created.
        int amount = 40;

        // -- layers
        // if the value is over 2, we can make layered cuboids. (set as default to 2)
        int layers = 1;
    };

int main() {
    if (!glfwInit())
        return -1;

    GLFWwindow* window = glfwCreateWindow(1600, 1200, "CGAL Viewer", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);

    // --- 設定 ---
    settings config;
    config.shape = 1; // 1-reflex (prism) ベースの生成モード
    config.mode = 1;  // 0: Intersectionチェックなし (shape = 1 にすることで 2-reflex が生成可能)
                      // 1: Intersectionチェックあり (shape = 1 にすることで 1-reflex が生成可能)
                      // ※ 2-reflexを作りたい場合はここを0にします！
    
    // --- 2-Reflex 生成用変数 ---
    double init_height = get_rand(100, 500, 800, mt_engine); // 固定高さ
    double min_x = 10000, max_x = -10000;
    double min_z = 10000, max_z = -10000;

    Nef_polyhedron cubes;

    // --- 最初の1つ ---
    // shape=1 なので、高さは固定(init_height)
    // 2-Reflex (浮いた島) にするために Y = init_height に設定
    // 1-Reflex にしたい場合は Y = 0.0 に変更してください
    double start_y = (config.mode == 0) ? init_height : 0.0;
    
    double len = get_rand(100, 500, 800, mt_engine);
    double dep = get_rand(100, 500, 800, mt_engine);
    Nef_polyhedron init_cube = make_cube(Point(0.0, start_y, 0.0), len, init_height, dep);
    
    cubes += init_cube;
    
    // 範囲記録 (2-reflexの土台用)
    min_x = -len/2; max_x = len/2; min_z = -dep/2; max_z = dep/2;

    std::cout << "Generating " << (config.mode == 0 ? "2-Reflex" : "1-Reflex") << " shape..." << std::endl;

    // --- ループ生成 ---
    int rate = 0; // 範囲を広げる用
    for(int i = 0; i < config.amount; i++){
        double x = get_rand(100, 500+rate, 1000+rate, mt_engine) * xyz(mt_engine);
        double z = get_rand(100, 500+rate, 1000+rate, mt_engine) * xyz(mt_engine);
        
        len = get_rand(100, 300+rate*0.5, 800+rate*0.5, mt_engine);
        dep = get_rand(100, 300+rate*0.5, 800+rate*0.5, mt_engine);

        // 範囲更新
        if (x - len/2 < min_x) min_x = x - len/2;
        if (x + len/2 > max_x) max_x = x + len/2;
        if (z - dep/2 < min_z) min_z = z - dep/2;
        if (z + dep/2 > max_z) max_z = z + dep/2;

        // 生成: Y座標は start_y (0 または init_height)
        Nef_polyhedron new_shape = make_cube(Point(x, start_y, z), len, init_height, dep);

        if (config.mode == 1){ 
            // 1-Reflex: !Intersectionチェックあり
            if (!cubes.intersection(new_shape).is_empty()) {
                cubes += new_shape;
                rate += 50;
            } else {
                std::cout << "Skipped (No intersection)" << std::endl;
            }
        }
        else {
            // 2-Reflex: Intersectionチェックあり
            if (cubes.intersection(new_shape).is_empty()) {
                cubes += new_shape;
                rate += 50;
            }
        }
    }

    // --- Bottom Cuboid (土台) の追加 (2-Reflexの場合のみ) ---
    if (config.shape == 1 && config.mode == 0) {
        std::cout << "Adding bottom plate for 2-Reflex..." << std::endl;
        // 全体をカバーする土台
        double base_width = (max_x - min_x); 
        double base_depth = (max_z - min_z);
        double center_x = (max_x + min_x) / 2.0;
        double center_z = (max_z + min_z) / 2.0;

        // Y=0 の位置に、上の図形と接するように配置
        Nef_polyhedron bottom_cube = make_cube(Point(center_x, 0.0, center_z), base_width, init_height, base_depth);
        cubes += bottom_cube;
    }

    Surface_mesh final_mesh;
    CGAL::convert_nef_polyhedron_to_polygon_mesh(cubes, final_mesh);

    tag_original_edges(final_mesh);

    // Added this 2 lines, for making triangulation and drawing certainly.
    assign_random_face_colors(final_mesh);
    CGAL::Polygon_mesh_processing::triangulate_faces(final_mesh);

    // いったん削除
    triangulate_faces_with_colors(final_mesh);

    std::cout << "\nClassifying edges into 24 patterns..." << std::endl;
    auto edge_pattern_map = classify_edges_into_patterns(final_mesh);
    std::cout << "Edge classification complete." << std::endl;

    // --- 2. 辺の総数と、各パターンの出現回数を集計 ---
    size_t real_edge_count = 0;
    size_t convex_edge_count = 0;
    size_t reflex_edge_count = 0;
    std::map<CornerPattern, int> pattern_counts;

    for (const auto& pair : edge_pattern_map) {
        CornerPattern pattern = pair.second;
        pattern_counts[pattern]++;
        
        if (pattern != CornerPattern::UNDEFINED) {
            real_edge_count++;
            if (static_cast<int>(pattern) >= static_cast<int>(CornerPattern::REFLEX_PX_PY)) {
                reflex_edge_count++;
            } else {
                convex_edge_count++;
            }
        }
    }

    std::cout << "\n--- Total Edge Counts ---" << std::endl;
    std::cout << "Total Real Edges (m): " << real_edge_count << std::endl;
    std::cout << "  - Convex Edges:     " << convex_edge_count << std::endl;
    std::cout << "  - Reflex Edges:     " << reflex_edge_count << std::endl;
    std::cout << "---------------------------\n" << std::endl;

    // --- 3. 12ケースの m/6 検証を実行 ---
    const auto target_sets = define_target_pattern_sets();
    size_t m = real_edge_count;
    double m_over_6 = static_cast<double>(m) / 6.0;

    std::cout << "--- m/6 Theorem Verification for 12 Sets ---" << std::endl;
    std::cout << "m = " << m << ", m/6 = " << m_over_6 << std::endl;
    std::cout << "--------------------------------------------" << std::endl;

    int min_guards = -1;
    int best_set_index = -1;
    for (int i = 0; i < target_sets.size(); ++i) {
        const auto& current_set = target_sets[i];
        if (current_set.size() != 3) continue; // 念のため
        
        auto it = current_set.begin();
        CornerPattern p1 = *it; ++it;
        CornerPattern p2 = *it; ++it;
        CornerPattern p3 = *it;

        int count1 = pattern_counts.count(p1) ? pattern_counts.at(p1) : 0;
        int count2 = pattern_counts.count(p2) ? pattern_counts.at(p2) : 0;
        int count3 = pattern_counts.count(p3) ? pattern_counts.at(p3) : 0;
        
        int guards_G = count1 + count2 + count3;

        if (best_set_index == -1 || guards_G < min_guards) {
            min_guards = guards_G;
            best_set_index = i;
        }
        
        if (i == 0)  std::cout << "--- XY Plane ---" << std::endl;
        if (i == 4)  std::cout << "--- YZ Plane ---" << std::endl;
        if (i == 8)  std::cout << "--- ZX Plane ---" << std::endl;
        
        std::cout << "Set " << i + 1 << ": G = " << guards_G;
        if (guards_G <= m_over_6) {
            std::cout << " (SUCCESS: G <= m/6) (" << count1 << ", " << count2 << ", " << count3 << ")"<< std::endl;
        } else {
            std::cout << " (SUCCESS: G <= m/6) (" << count1 << ", " << count2 << ", " << count3 << ")"<< std::endl;
        }
    }
    std::cout << "--------------------------------------------\n" << std::endl;
    // ステップ4: 最も効率の良かったセットの辺をハイライトする
    if (best_set_index != -1) {
        std::cout << "Highlighting edges for the best set (Set #" << best_set_index + 1 << ")..." << std::endl;
        
        // best_set_indexを使って、ターゲットセットの中から最良のセットを取得
        const std::set<CornerPattern>& best_set_patterns = target_sets[best_set_index];
        
        // タグ付け関数を呼び出し
        tag_highlight_edges(final_mesh, edge_pattern_map, best_set_patterns);
    }
    
    while (!glfwWindowShouldClose(window)) {
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        float aspect = (float)width / (float)height;

        glViewport(0, 0, width, height);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        set_perspective(45.0f, aspect, 0.1f, 10.0f);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glTranslatef(0.0f, 0.0f, -3.0f);
        glRotatef(rotation_x, 1.0f, 0.0f, 0.0f);
        glRotatef(rotation_y, 0.0f, 1.0f, 0.0f);

        draw_mesh(final_mesh);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}

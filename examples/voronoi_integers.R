#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <boost/polygon/voronoi.hpp>
#include <cmath>
#include <limits>
#include <algorithm>
#include <set>

namespace bp = boost::polygon;
using namespace std;

// Define a simple point structure
struct Point {
  double x, y;
};

// Specialization for Boost Polygon
namespace boost {
  namespace polygon {
    template <>
      struct geometry_concept<Point> {
        typedef point_concept type;
      };
    template <>
      struct point_traits<Point> {
        typedef double coordinate_type;
        static inline coordinate_type get(const Point& p, orientation_2d orient) {
          return (orient == HORIZONTAL) ? p.x : p.y;
        }
      };
  }  // namespace polygon
}  // namespace boost

// [[Rcpp::export]]
Rcpp::List generate_voronoi(Rcpp::NumericMatrix points,
                            double x_min, double x_max, double y_min, double y_max) {
  typedef bp::voronoi_diagram<double> VoronoiDiagram;
  VoronoiDiagram vd;

  std::vector<Point> point_vec;

  double x_smallest_value = std::numeric_limits<double>::infinity();
  double x_largest_value  = -std::numeric_limits<double>::infinity();
  double y_smallest_value = std::numeric_limits<double>::infinity();
  double y_largest_value  = -std::numeric_limits<double>::infinity();

  for (int i = 0; i < points.nrow(); i++) {
    double x = points(i, 0);
    double y = points(i, 1);

    // Track smallest/largest values
    if (x < x_smallest_value) x_smallest_value = x;
    if (x > x_largest_value)  x_largest_value = x;
    if (y < y_smallest_value) y_smallest_value = y;
    if (y > y_largest_value)  y_largest_value = y;

    // Adapt x/y limits if point is outside
    if (x < x_min) {
      Rcpp::Rcout << "Input point x < x_min. Updating x_min <- " << x << std::endl;
      x_min = x;
    }
    if (x > x_max) {
      Rcpp::Rcout << "Input point x > x_max. Updating x_max <- " << x << std::endl;
      x_max = x;
    }
    if (y < y_min) {
      Rcpp::Rcout << "Input point y < y_min. Updating y_min <- " << y << std::endl;
      y_min = y;
    }
    if (y > y_max) {
      Rcpp::Rcout << "Input point y > y_max. Updating y_max <- " << y << std::endl;
      y_max = y;
    }

    point_vec.push_back({x, y});
  }



  bp::construct_voronoi(point_vec.begin(), point_vec.end(), &vd);

  Rcpp::List finite_edges;
  Rcpp::List infinite_edges;
  Rcpp::List site3_list;
  std::vector<std::pair<double, double>> inf_edges;
  std::vector<std::vector<std::vector<double>>> cell_polygons(point_vec.size());

  const double tol = 1e-8;
  std::map<const bp::voronoi_edge<double>*, std::pair<double, double>> infinite_endpoints;

  for (const auto& edge : vd.edges()) {
    if (edge.vertex0() && edge.vertex1()) {
      double x1 = edge.vertex0()->x();

      double y1 = edge.vertex0()->y();
      double x2 = edge.vertex1()->x();
      double y2 = edge.vertex1()->y();
      finite_edges.push_back(Rcpp::NumericVector::create(x1, y1, x2, y2));

      // Update bounding box based on finite edge vertices
      if (x1 < x_min){
        Rcpp::Rcout << "Finite edge vertices were beyond the bounding box limits." << std::endl;
        Rcpp::Rcout << "Updating x_min <- " << x1 << std::endl;
        x_min = x1;
      }
      if (x2 < x_min){
        Rcpp::Rcout << "Finite edge vertices were beyond the bounding box limits." << std::endl;
        Rcpp::Rcout << "Updating x_min <- " << x2 << std::endl;
        x_min = x2;
      }
      if (x1 > x_max){
        Rcpp::Rcout << "Finite edge vertices were beyond the bounding box limits." << std::endl;
        Rcpp::Rcout << "Updating x_max <- " << x1 << std::endl;
        x_max = x1;
      }
      if (x2 > x_max){
        Rcpp::Rcout << "Finite edge vertices were beyond the bounding box limits." << std::endl;
        Rcpp::Rcout << "Updating x_max <- " << x2 << std::endl;
        x_max = x2;
      }
      if (y1 < y_min){
        Rcpp::Rcout << "Finite edge vertices were beyond the bounding box limits." << std::endl;
        Rcpp::Rcout << "Updating y_min <- " << y1 << std::endl;
        y_min = y1;
      }
      if (y2 < y_min){
        Rcpp::Rcout << "Finite edge vertices were beyond the bounding box limits." << std::endl;
        Rcpp::Rcout << "Updating y_min <- " << y2 << std::endl;
        y_min = y2;
      }
      if (y1 > y_max){
        Rcpp::Rcout << "Finite edge vertices were beyond the bounding box limits." << std::endl;
        Rcpp::Rcout << "Updating y_max <- " << y1 << std::endl;
        y_max = y1;
      }
      if (y2 > y_max) {
        Rcpp::Rcout << "Finite edge vertices were beyond the bounding box limits." << std::endl;
        Rcpp::Rcout << "Updating y_max <- " << y2 << std::endl;
        y_max = y2;
      }

    }
  }
  for (const auto& edge : vd.edges()) {
    if (!(edge.vertex0() && edge.vertex1())) {
      const bp::voronoi_cell<double>* c1 = edge.cell();
      const bp::voronoi_cell<double>* c2 = edge.twin()->cell();
      if (!c1 || !c2) continue;

      int idx1 = c1->source_index();
      int idx2 = c2->source_index();
      Point site1 = point_vec[idx1];
      Point site2 = point_vec[idx2];

      double vx, vy;
      if (edge.vertex0()) {
        vx = edge.vertex0()->x();
        vy = edge.vertex0()->y();
      } else {
        vx = edge.vertex1()->x();
        vy = edge.vertex1()->y();
      }

      double r_sq = (site1.x - vx) * (site1.x - vx) + (site1.y - vy) * (site1.y - vy);

      int site3_index = -1;
      double x3 = std::numeric_limits<double>::quiet_NaN();
      double y3 = std::numeric_limits<double>::quiet_NaN();
      for (std::size_t i = 0; i < point_vec.size(); i++) {
        if (i == static_cast<std::size_t>(idx1) || i == static_cast<std::size_t>(idx2)) continue;
        Point p = point_vec[i];
        double d_sq = (p.x - vx) * (p.x - vx) + (p.y - vy) * (p.y - vy);
        if (std::abs(d_sq - r_sq) < tol) {
          site3_index = static_cast<int>(i);
          x3 = p.x;
          y3 = p.y;
          break;
        }
      }

      double dx = site2.y - site1.y;
      double dy = site1.x - site2.x;
      double len = sqrt(dx * dx + dy * dy);
      dx /= len;
      dy /= len;

      int direction_sign = 1;
      if (site3_index >= 0) {
        double dot_to_site3 = (x3 - vx) * dx + (y3 - vy) * dy;
        if (dot_to_site3 > 0) {
          direction_sign = -1;
        }
      }

      double t_x = std::numeric_limits<double>::infinity();
      double t_y = std::numeric_limits<double>::infinity();

      if (dx != 0) {
        t_x = (direction_sign * dx > 0)
        ? (x_max - vx) / (direction_sign * dx)
        : (x_min - vx) / (direction_sign * dx);
      }
      if (dy != 0) {
        t_y = (direction_sign * dy > 0)
        ? (y_max - vy) / (direction_sign * dy)
        : (y_min - vy) / (direction_sign * dy);
      }

      double t = std::min(t_x, t_y);
      double x2 = vx + direction_sign * dx * t;
      double y2 = vy + direction_sign * dy * t;

      infinite_edges.push_back(Rcpp::NumericVector::create(vx, vy, x2, y2));
      inf_edges.push_back({x2, y2});
      infinite_endpoints[&edge] = std::make_pair(x2, y2);

      if (site3_index >= 0) {
        site3_list.push_back(Rcpp::NumericVector::create(x3, y3));
      } else {
        site3_list.push_back(R_NilValue);
      }
    }
  }

  // Boundary edges
  std::vector<Rcpp::NumericVector> boundary_edges;
  auto filter_and_sort = [&](bool vertical, double val) {
    std::vector<double> coords;
    for (auto& p : inf_edges) {
      if (vertical) {
        if (std::abs(p.first - val) < tol) coords.push_back(p.second);
      } else {
        if (std::abs(p.second - val) < tol) coords.push_back(p.first);
      }
    }
    std::sort(coords.begin(), coords.end());
    return coords;
  };

  auto y_coords = filter_and_sort(true, x_min);
  double prev_y = y_min;
  for (double y : y_coords) {
    boundary_edges.push_back(Rcpp::NumericVector::create(x_min, prev_y, x_min, y));
    prev_y = y;
  }
  boundary_edges.push_back(Rcpp::NumericVector::create(x_min, prev_y, x_min, y_max));

  auto x_coords = filter_and_sort(false, y_max);
  double prev_x = x_min;
  for (double x : x_coords) {
    boundary_edges.push_back(Rcpp::NumericVector::create(prev_x, y_max, x, y_max));
    prev_x = x;
  }
  boundary_edges.push_back(Rcpp::NumericVector::create(prev_x, y_max, x_max, y_max));

  y_coords = filter_and_sort(true, x_max);
  std::reverse(y_coords.begin(), y_coords.end());
  prev_y = y_max;
  for (double y : y_coords) {
    boundary_edges.push_back(Rcpp::NumericVector::create(x_max, prev_y, x_max, y));
    prev_y = y;
  }
  boundary_edges.push_back(Rcpp::NumericVector::create(x_max, prev_y, x_max, y_min));

  x_coords = filter_and_sort(false, y_min);
  std::reverse(x_coords.begin(), x_coords.end());
  prev_x = x_max;
  for (double x : x_coords) {
    boundary_edges.push_back(Rcpp::NumericVector::create(prev_x, y_min, x, y_min));
    prev_x = x;
  }
  boundary_edges.push_back(Rcpp::NumericVector::create(prev_x, y_min, x_min, y_min));

  // Build polygons
  for (const auto& cell : vd.cells()) {
    std::vector<std::vector<double>> polygon;
    const auto* edge = cell.incident_edge();
    const auto* start = edge;

    do {
      if (edge->is_primary()) {
        if (edge->vertex0()) {
          polygon.push_back({edge->vertex0()->x(), edge->vertex0()->y()});
        } else if (infinite_endpoints.count(edge)) {
          polygon.push_back({infinite_endpoints[edge].first, infinite_endpoints[edge].second});
        }

        if (edge->vertex1() == nullptr && infinite_endpoints.count(edge)) {
          polygon.push_back({infinite_endpoints[edge].first, infinite_endpoints[edge].second});
        }
      }
      edge = edge->next();
    } while (edge != start);
    cell_polygons[cell.source_index()] = polygon;
  }

  // Site coordinates
  Rcpp::List polygon_sites;
  for (const auto& pt : point_vec) {
    polygon_sites.push_back(Rcpp::NumericVector::create(pt.x, pt.y));
  }

  // Add the boundry vertices
  for (size_t i = 0; i < cell_polygons.size(); i++) {

    int on_x_min_limit = 0;
    int on_x_max_limit = 0;
    int on_y_min_limit = 0;
    int on_y_max_limit = 0;

    for (const auto& vertex : cell_polygons[i]) {
      double x = vertex[0];
      double y = vertex[1];

      // Check if vertex is on the bounding box
      if (std::abs(x - x_min) < tol){
        on_x_min_limit++;
      }
      if (std::abs(x - x_max) < tol){
        on_x_max_limit++;
      }
      if (std::abs(y - y_min) < tol){
        on_y_min_limit++;
      }
      if (std::abs(y - y_max) < tol){
        on_y_max_limit++;
      }
    }
    if(on_x_min_limit >= 1 && on_y_min_limit >= 1){
      cell_polygons[i].push_back({x_min, y_min});
    }
    if(on_x_min_limit >= 1 && on_y_max_limit >= 1){
      cell_polygons[i].push_back({x_min, y_max});
    }
    if(on_x_max_limit >= 1 && on_y_min_limit >= 1){
      cell_polygons[i].push_back({x_max, y_min});
    }
    if(on_x_max_limit >= 1 && on_y_max_limit >= 1){
      cell_polygons[i].push_back({x_max, y_max});
    }
    if(on_x_min_limit >= 1 && on_x_max_limit >= 1){
      Rcpp::NumericVector site = polygon_sites[i];
      double y_value = site[1];

      if (y_value >= y_largest_value){
        cell_polygons[i].push_back({x_min, y_max});
        cell_polygons[i].push_back({x_max, y_max});
      } else if(y_value <= y_smallest_value){
        cell_polygons[i].push_back({x_min, y_min});
        cell_polygons[i].push_back({x_max, y_min});
      }
    }

    if(on_y_min_limit >= 1 && on_y_max_limit >= 1){
      Rcpp::NumericVector site = polygon_sites[i];
      double x_value = site[0];

      if (x_value >= x_largest_value){
        cell_polygons[i].push_back({x_max, y_min});
        cell_polygons[i].push_back({x_max, y_max});
      } else if (x_value <= x_smallest_value){
        cell_polygons[i].push_back({x_min, y_min});
        cell_polygons[i].push_back({x_min, y_max});
      }
    }
  }

  // Sort polygon vertices counter-clockwise
  for (auto& polygon : cell_polygons) {
    int n = polygon.size();
    if (n > 2) {
      // Compute centroid
      double cx = 0.0, cy = 0.0;
      for (const auto& v : polygon) {
        cx += v[0];
        cy += v[1];
      }
      cx /= n;
      cy /= n;

      // Sort vertices by angle from centroid
      std::sort(polygon.begin(), polygon.end(), [&](const std::vector<double>& a, const std::vector<double>& b) {
        double angle_a = std::atan2(a[1] - cy, a[0] - cx);
        double angle_b = std::atan2(b[1] - cy, b[0] - cx);
        return angle_a < angle_b;
      });
    }
  }


  // Polygon areas
  std::vector<double> polygon_areas;
  for (const auto& polygon : cell_polygons) {
    double area = 0.0;
    int n = polygon.size();
    if (n > 2) {
      for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        area += polygon[i][0] * polygon[j][1] - polygon[j][0] * polygon[i][1];
      }
      area = std::abs(area) * 0.5;
    }
    polygon_areas.push_back(area);
  }

  // Collect all unique vertices
  auto make_vertex_key = [&](double x, double y) {
    int scale = 1e8;
    return std::make_pair((int)(x * scale), (int)(y * scale));
  };

  std::set<std::pair<int, int>> vertex_keys;
  std::vector<Rcpp::NumericVector> all_vertices;

  auto add_vertex = [&](double x, double y) {
    auto key = make_vertex_key(x, y);
    if (vertex_keys.insert(key).second) {
      all_vertices.push_back(Rcpp::NumericVector::create(x, y));
    }
  };

  for (auto& edge : finite_edges) {
    Rcpp::NumericVector e = edge;
    add_vertex(e[0], e[1]);
    add_vertex(e[2], e[3]);
  }

  for (auto& edge : infinite_edges) {
    Rcpp::NumericVector e = edge;
    add_vertex(e[0], e[1]);
    add_vertex(e[2], e[3]);
  }

  for (auto& edge : boundary_edges) {
    Rcpp::NumericVector e = edge;
    add_vertex(e[0], e[1]);
    add_vertex(e[2], e[3]);
  }



  // Return the indices of polygons touching at least two boundary points
  return Rcpp::List::create(
    Rcpp::Named("site_points") = polygon_sites,
    Rcpp::Named("finite_edges") = finite_edges,
    Rcpp::Named("infinite_edges") = infinite_edges,
    Rcpp::Named("site3") = site3_list,
    Rcpp::Named("boundary_edges") = boundary_edges,
    Rcpp::Named("polygons") = cell_polygons,
    Rcpp::Named("polygon_areas") = polygon_areas,
    Rcpp::Named("vertices") = all_vertices,
    Rcpp::Named("x_min") = x_min,
    Rcpp::Named("x_max") = x_max,
    Rcpp::Named("y_min") = y_min,
    Rcpp::Named("y_max") = y_max
  );
}

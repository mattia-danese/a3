#ifndef SPHERE_H
#define SPHERE_H

#include <vector>
#include <map>
#include <set>
#include <cmath>
#include "Shape.h"

using namespace std;

class Sphere : public Shape {
public:
	Sphere() {};
	~Sphere() {};

	OBJ_TYPE getType() {
		return SHAPE_SPHERE;
	}

	void draw() {
		float m_segmentsX_angle = 2 * PI / m_segmentsX;
		float m_segmentsY_angle = PI / m_segmentsY;

		normals_by_vertex.clear();
		map<vector<float>, vector<vector<float>>> side_points = {};

		for (int x_index = 0; x_index < m_segmentsX; x_index++) {
			for (int y_index = 0; y_index < m_segmentsY; y_index++) {
				vector<float> each_side_normal = getNormal(0.5f * cos(PI / 2.0f - m_segmentsY_angle * y_index) * sin(m_segmentsX_angle * x_index), 0.5f * sin(PI / 2.0f - m_segmentsY_angle * y_index), 0.5f * cos(PI / 2.0f - m_segmentsY_angle * y_index) * cos(m_segmentsX_angle * x_index), 0.5f * cos(PI / 2.0f - m_segmentsY_angle * (y_index + 1)) * sin(m_segmentsX_angle * x_index), 0.5f * sin(PI / 2.0f - m_segmentsY_angle * (y_index + 1)), 0.5f * cos(PI / 2.0f - m_segmentsY_angle * (y_index + 1)) * cos(m_segmentsX_angle * x_index), 0.5f * cos(PI / 2.0f - m_segmentsY_angle * y_index) * sin(m_segmentsX_angle * (x_index + 1)), 0.5f * sin(PI / 2.0f - m_segmentsY_angle * y_index), 0.5f * cos(PI / 2.0f - m_segmentsY_angle * y_index) * cos(m_segmentsX_angle * (x_index + 1)));
				each_side_normal = getNormalizedVector(each_side_normal[0], each_side_normal[1], each_side_normal[2]);
				vector<vector<float>> each_side_points = {};

				each_side_points.push_back({ 0.5f * cos(PI / 2.0f - m_segmentsY_angle * y_index) * sin(m_segmentsX_angle * x_index), 0.5f * sin(PI / 2.0f - m_segmentsY_angle * y_index), 0.5f * cos(PI / 2.0f - m_segmentsY_angle * y_index) * cos(m_segmentsX_angle * x_index) });
				each_side_points.push_back({ 0.5f * cos(PI / 2.0f - m_segmentsY_angle * (y_index + 1)) * sin(m_segmentsX_angle * x_index), 0.5f * sin(PI / 2.0f - m_segmentsY_angle * (y_index + 1)), 0.5f * cos(PI / 2.0f - m_segmentsY_angle * (y_index + 1)) * cos(m_segmentsX_angle * x_index) });
				each_side_points.push_back({ 0.5f * cos(PI / 2.0f - m_segmentsY_angle * y_index) * sin(m_segmentsX_angle * (x_index + 1)), 0.5f * sin(PI / 2.0f - m_segmentsY_angle * y_index), 0.5f * cos(PI / 2.0f - m_segmentsY_angle * y_index) * cos(m_segmentsX_angle * (x_index + 1)) });

				normals_by_vertex[{ 0.5f * cos(PI / 2.0f - m_segmentsY_angle * y_index) * sin(m_segmentsX_angle * x_index), 0.5f * sin(PI / 2.0f - m_segmentsY_angle * y_index), 0.5f * cos(PI / 2.0f - m_segmentsY_angle * y_index) * cos(m_segmentsX_angle * x_index) }].insert(getNormalizedVector(0.5f * cos(PI / 2.0f - m_segmentsY_angle * y_index) * sin(m_segmentsX_angle * x_index), 0.5f * sin(PI / 2.0f - m_segmentsY_angle * y_index), 0.5f * cos(PI / 2.0f - m_segmentsY_angle * y_index) * cos(m_segmentsX_angle * x_index)));

				side_points[each_side_normal].insert(side_points[each_side_normal].end(), each_side_points.begin(), each_side_points.end());

				each_side_normal = getNormal(0.5f * cos(PI / 2.0f - m_segmentsY_angle * (y_index + 1)) * sin(m_segmentsX_angle * x_index), 0.5f * sin(PI / 2.0f - m_segmentsY_angle * (y_index + 1)), 0.5f * cos(PI / 2.0f - m_segmentsY_angle * (y_index + 1)) * cos(m_segmentsX_angle * x_index), 0.5f * cos(PI / 2.0f - m_segmentsY_angle * (y_index + 1)) * sin(m_segmentsX_angle * (x_index + 1)), 0.5f * sin(PI / 2.0f - m_segmentsY_angle * (y_index + 1)), 0.5f * cos(PI / 2.0f - m_segmentsY_angle * (y_index + 1)) * cos(m_segmentsX_angle * (x_index + 1)), 0.5f * cos(PI / 2.0f - m_segmentsY_angle * y_index) * sin(m_segmentsX_angle * (x_index + 1)), 0.5f * sin(PI / 2.0f - m_segmentsY_angle * y_index), 0.5f * cos(PI / 2.0f - m_segmentsY_angle * y_index) * cos(m_segmentsX_angle * (x_index + 1)));
				each_side_normal = getNormalizedVector(each_side_normal[0], each_side_normal[1], each_side_normal[2]);
				each_side_points.clear();
				each_side_points.push_back({ 0.5f * cos(PI / 2.0f - m_segmentsY_angle * (y_index + 1)) * sin(m_segmentsX_angle * x_index), 0.5f * sin(PI / 2.0f - m_segmentsY_angle * (y_index + 1)), 0.5f * cos(PI / 2.0f - m_segmentsY_angle * (y_index + 1)) * cos(m_segmentsX_angle * x_index) });
				each_side_points.push_back({ 0.5f * cos(PI / 2.0f - m_segmentsY_angle * (y_index + 1)) * sin(m_segmentsX_angle * (x_index + 1)), 0.5f * sin(PI / 2.0f - m_segmentsY_angle * (y_index + 1)), 0.5f * cos(PI / 2.0f - m_segmentsY_angle * (y_index + 1)) * cos(m_segmentsX_angle * (x_index + 1)) });
				each_side_points.push_back({ 0.5f * cos(PI / 2.0f - m_segmentsY_angle * y_index) * sin(m_segmentsX_angle * (x_index + 1)), 0.5f * sin(PI / 2.0f - m_segmentsY_angle * y_index), 0.5f * cos(PI / 2.0f - m_segmentsY_angle * y_index) * cos(m_segmentsX_angle * (x_index + 1)) });

				normals_by_vertex[{ 0.5f * cos(PI / 2.0f - m_segmentsY_angle * (y_index + 1)) * sin(m_segmentsX_angle * x_index), 0.5f * sin(PI / 2.0f - m_segmentsY_angle * (y_index + 1)), 0.5f * cos(PI / 2.0f - m_segmentsY_angle * (y_index + 1)) * cos(m_segmentsX_angle * x_index) }].insert(getNormalizedVector(0.5f * cos(PI / 2.0f - m_segmentsY_angle * (y_index + 1)) * sin(m_segmentsX_angle * x_index), 0.5f * sin(PI / 2.0f - m_segmentsY_angle * (y_index + 1)), 0.5f * cos(PI / 2.0f - m_segmentsY_angle * (y_index + 1)) * cos(m_segmentsX_angle * x_index)));

				side_points[each_side_normal].insert(side_points[each_side_normal].end(), each_side_points.begin(), each_side_points.end());
			}
		}

		glBegin(GL_TRIANGLES);
		for (map<vector<float>, vector<vector<float>>>::iterator iter = side_points.begin(); iter != side_points.end(); iter++) {
			for (int index = 0; index < iter->second.size(); index++) {
					normalizeNormal(iter->second[index][0], iter->second[index][1], iter->second[index][2]);
					glVertex3f(iter->second[index][0], iter->second[index][1], iter->second[index][2]);
			}
		}
		glEnd();
	};

	void drawNormal() {
		float length_factor = 0.1f;

		glBegin(GL_LINES);
		glColor3f(1.0, 0.0, 0.0);
		for (map<vector<float>, set<vector<float>>>::iterator iter = normals_by_vertex.begin(); iter != normals_by_vertex.end(); iter++) {
			for (auto it = iter->second.cbegin(); it != iter->second.cend(); it++) {
				glVertex3f(iter->first[0], iter->first[1], iter->first[2]);
				glVertex3f(iter->first[0] + (*it)[0] * length_factor, iter->first[1] + (*it)[1] * length_factor, iter->first[2] + (*it)[2] * length_factor);
			}
		}
		glEnd();
	};

private:
	map<vector<float>, set<vector<float>>> normals_by_vertex;
};

#endif
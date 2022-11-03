#ifndef CUBE_H
#define CUBE_H

#include <vector>
#include <map>
#include <set>
#include "Shape.h"

using namespace std;

class Cube : public Shape {
public:
	Cube() {};
	~Cube() {};

	OBJ_TYPE getType() {
		return SHAPE_CUBE;
	}

	void draw() {
		float m_segmentsX_length = 1.0f / m_segmentsX;
		float m_segmentsY_length = 1.0f / m_segmentsY;

		normals_by_vertex.clear();
		vector<vector<float>> front_points = {};

		for (int x_index = 0; x_index < m_segmentsX; x_index++) {
			for (int y_index = 0; y_index < m_segmentsY; y_index++) {
				front_points.push_back({ -0.5f + m_segmentsX_length * x_index, -0.5f + m_segmentsY_length * y_index, 0.5f, 1.0f });
				front_points.push_back({ -0.5f + m_segmentsX_length * (x_index + 1), -0.5f + m_segmentsY_length * y_index, 0.5f, 1.0f });
				front_points.push_back({ -0.5f + m_segmentsX_length * (x_index + 1), -0.5f + m_segmentsY_length * (y_index + 1), 0.5f, 1.0f });

				front_points.push_back({ -0.5f + m_segmentsX_length * x_index, -0.5f + m_segmentsY_length * y_index, 0.5f, 1.0f });
				front_points.push_back({ -0.5f + m_segmentsX_length * (x_index + 1), -0.5f + m_segmentsY_length * (y_index + 1), 0.5f, 1.0f }); 
				front_points.push_back({ -0.5f + m_segmentsX_length * x_index, -0.5f + m_segmentsY_length * (y_index + 1), 0.5f, 1.0f });

				normals_by_vertex[{ -0.5f + m_segmentsX_length * x_index, -0.5f + m_segmentsY_length * y_index, 0.5f, 1.0f }].insert({ 0.0f, 0.0f, 1.0f });
				normals_by_vertex[{ -0.5f + m_segmentsX_length * (x_index + 1), -0.5f + m_segmentsY_length * y_index, 0.5f, 1.0f }].insert({ 0.0f, 0.0f, 1.0f });
				normals_by_vertex[{ -0.5f + m_segmentsX_length * (x_index + 1), -0.5f + m_segmentsY_length * (y_index + 1), 0.5f, 1.0f }].insert({ 0.0f, 0.0f, 1.0f });

				normals_by_vertex[{ -0.5f + m_segmentsX_length * x_index, -0.5f + m_segmentsY_length * y_index, 0.5f, 1.0f }].insert({ 0.0f, 0.0f, 1.0f });
				normals_by_vertex[{ -0.5f + m_segmentsX_length * (x_index + 1), -0.5f + m_segmentsY_length * (y_index + 1), 0.5f, 1.0f }].insert({ 0.0f, 0.0f, 1.0f });
				normals_by_vertex[{ -0.5f + m_segmentsX_length * x_index, -0.5f + m_segmentsY_length * (y_index + 1), 0.5f, 1.0f }].insert({ 0.0f, 0.0f, 1.0f });
			}
		}

		vector<vector<float>> back_points = {};

		for (int x_index = 0; x_index < m_segmentsX; x_index++) {
			for (int y_index = 0; y_index < m_segmentsY; y_index++) {
				back_points.push_back({ -0.5f + m_segmentsX_length * (x_index + 1), -0.5f + m_segmentsY_length * y_index, -0.5f, 1.0f });
				back_points.push_back({ -0.5f + m_segmentsX_length * x_index, -0.5f + m_segmentsY_length * y_index, -0.5f, 1.0f });
				back_points.push_back({ -0.5f + m_segmentsX_length * x_index, -0.5f + m_segmentsY_length * (y_index + 1), -0.5f, 1.0f });

				back_points.push_back({ -0.5f + m_segmentsX_length * (x_index + 1), -0.5f + m_segmentsY_length * y_index, -0.5f, 1.0f });
				back_points.push_back({ -0.5f + m_segmentsX_length * x_index, -0.5f + m_segmentsY_length * (y_index + 1), -0.5f, 1.0f });
				back_points.push_back({ -0.5f + m_segmentsX_length * (x_index + 1), -0.5f + m_segmentsY_length * (y_index + 1), -0.5f, 1.0f });

				normals_by_vertex[{ -0.5f + m_segmentsX_length * (x_index + 1), -0.5f + m_segmentsY_length * y_index, -0.5f, 1.0f }].insert({ 0.0f, 0.0f, -1.0f });
				normals_by_vertex[{ -0.5f + m_segmentsX_length * x_index, -0.5f + m_segmentsY_length * y_index, -0.5f, 1.0f }].insert({ 0.0f, 0.0f, -1.0f });
				normals_by_vertex[{ -0.5f + m_segmentsX_length * x_index, -0.5f + m_segmentsY_length * (y_index + 1), -0.5f, 1.0f }].insert({ 0.0f, 0.0f, -1.0f });

				normals_by_vertex[{ -0.5f + m_segmentsX_length * (x_index + 1), -0.5f + m_segmentsY_length * y_index, -0.5f, 1.0f }].insert({ 0.0f, 0.0f, -1.0f });
				normals_by_vertex[{ -0.5f + m_segmentsX_length * x_index, -0.5f + m_segmentsY_length * (y_index + 1), -0.5f, 1.0f }].insert({ 0.0f, 0.0f, -1.0f });
				normals_by_vertex[{ -0.5f + m_segmentsX_length * (x_index + 1), -0.5f + m_segmentsY_length * (y_index + 1), -0.5f, 1.0f }].insert({ 0.0f, 0.0f, -1.0f });
			}
		}

		vector<vector<float>> top_points = {};

		for (int x_index = 0; x_index < m_segmentsX; x_index++) {
			for (int y_index = 0; y_index < m_segmentsY; y_index++) {
				top_points.push_back({ -0.5f + m_segmentsX_length * x_index, 0.5f, -0.5f + m_segmentsY_length * y_index, 1.0f });
				top_points.push_back({ -0.5f + m_segmentsX_length * x_index, 0.5f, -0.5f + m_segmentsY_length * (y_index + 1), 1.0f });
				top_points.push_back({ -0.5f + m_segmentsX_length * (x_index + 1), 0.5f, -0.5f + m_segmentsY_length * y_index, 1.0f });

				top_points.push_back({ -0.5f + m_segmentsX_length * (x_index + 1), 0.5f, -0.5f + m_segmentsY_length * y_index, 1.0f });
				top_points.push_back({ -0.5f + m_segmentsX_length * x_index, 0.5f, -0.5f + m_segmentsY_length * (y_index + 1), 1.0f });
				top_points.push_back({ -0.5f + m_segmentsX_length * (x_index + 1), 0.5f, -0.5f + m_segmentsY_length * (y_index + 1), 1.0f });

				normals_by_vertex[{ -0.5f + m_segmentsX_length * x_index, 0.5f, -0.5f + m_segmentsY_length * y_index, 1.0f }].insert({ 0.0f, 1.0f, 0.0f });
				normals_by_vertex[{ -0.5f + m_segmentsX_length * x_index, 0.5f, -0.5f + m_segmentsY_length * (y_index + 1), 1.0f }].insert({ 0.0f, 1.0f, 0.0f });
				normals_by_vertex[{ -0.5f + m_segmentsX_length * (x_index + 1), 0.5f, -0.5f + m_segmentsY_length * y_index, 1.0f }].insert({ 0.0f, 1.0f, 0.0f });

				normals_by_vertex[{ -0.5f + m_segmentsX_length * (x_index + 1), 0.5f, -0.5f + m_segmentsY_length * y_index, 1.0f }].insert({ 0.0f, 1.0f, 0.0f });
				normals_by_vertex[{ -0.5f + m_segmentsX_length * x_index, 0.5f, -0.5f + m_segmentsY_length * (y_index + 1), 1.0f }].insert({ 0.0f, 1.0f, 0.0f });
				normals_by_vertex[{ -0.5f + m_segmentsX_length * (x_index + 1), 0.5f, -0.5f + m_segmentsY_length * (y_index + 1), 1.0f }].insert({ 0.0f, 1.0f, 0.0f });
			}
		}

		vector<vector<float>> bottom_points = {};

		for (int x_index = 0; x_index < m_segmentsX; x_index++) {
			for (int y_index = 0; y_index < m_segmentsY; y_index++) {
				bottom_points.push_back({ -0.5f + m_segmentsX_length * x_index, -0.5f, -0.5f + m_segmentsY_length * y_index, 1.0f });
				bottom_points.push_back({ -0.5f + m_segmentsX_length * (x_index + 1), -0.5f, -0.5f + m_segmentsY_length * y_index, 1.0f });
				bottom_points.push_back({ -0.5f + m_segmentsX_length * (x_index + 1), -0.5f, -0.5f + m_segmentsY_length * (y_index + 1), 1.0f });

				bottom_points.push_back({ -0.5f + m_segmentsX_length * x_index, -0.5f, -0.5f + m_segmentsY_length * y_index, 1.0f });
				bottom_points.push_back({ -0.5f + m_segmentsX_length * (x_index + 1), -0.5f, -0.5f + m_segmentsY_length * (y_index + 1), 1.0f });
				bottom_points.push_back({ -0.5f + m_segmentsX_length * x_index, -0.5f, -0.5f + m_segmentsY_length * (y_index + 1), 1.0f });

				normals_by_vertex[{ -0.5f + m_segmentsX_length * x_index, -0.5f, -0.5f + m_segmentsY_length * y_index, 1.0f }].insert({ 0.0f, -1.0f, 0.0f });
				normals_by_vertex[{ -0.5f + m_segmentsX_length * (x_index + 1), -0.5f, -0.5f + m_segmentsY_length * y_index, 1.0f }].insert({ 0.0f, -1.0f, 0.0f });
				normals_by_vertex[{ -0.5f + m_segmentsX_length * (x_index + 1), -0.5f, -0.5f + m_segmentsY_length * (y_index + 1), 1.0f }].insert({ 0.0f, -1.0f, 0.0f });

				normals_by_vertex[{ -0.5f + m_segmentsX_length * x_index, -0.5f, -0.5f + m_segmentsY_length * y_index, 1.0f }].insert({ 0.0f, -1.0f, 0.0f });
				normals_by_vertex[{ -0.5f + m_segmentsX_length * (x_index + 1), -0.5f, -0.5f + m_segmentsY_length * (y_index + 1), 1.0f }].insert({ 0.0f, -1.0f, 0.0f });
				normals_by_vertex[{ -0.5f + m_segmentsX_length * x_index, -0.5f, -0.5f + m_segmentsY_length * (y_index + 1), 1.0f }].insert({ 0.0f, -1.0f, 0.0f });
			}
		}

		vector<vector<float>> left_points = {};

		for (int x_index = 0; x_index < m_segmentsX; x_index++) {
			for (int y_index = 0; y_index < m_segmentsY; y_index++) {
				left_points.push_back({ -0.5f, -0.5f + m_segmentsY_length * y_index, -0.5f + m_segmentsX_length * x_index, 1.0f });
				left_points.push_back({ -0.5f, -0.5f + m_segmentsY_length * y_index, -0.5f + m_segmentsX_length * (x_index + 1), 1.0f });
				left_points.push_back({ -0.5f, -0.5f + m_segmentsY_length * (y_index + 1), -0.5f + m_segmentsX_length * (x_index + 1), 1.0f });

				left_points.push_back({ -0.5f, -0.5f + m_segmentsY_length * y_index, -0.5f + m_segmentsX_length * x_index, 1.0f });
				left_points.push_back({ -0.5f, -0.5f + m_segmentsY_length * (y_index + 1), -0.5f + m_segmentsX_length * (x_index + 1), 1.0f });
				left_points.push_back({ -0.5f, -0.5f + m_segmentsY_length * (y_index + 1), -0.5f + m_segmentsX_length * x_index, 1.0f });

				normals_by_vertex[{ -0.5f, -0.5f + m_segmentsY_length * y_index, -0.5f + m_segmentsX_length * x_index, 1.0f }].insert({ -1.0f, 0.0f, 0.0f });
				normals_by_vertex[{ -0.5f, -0.5f + m_segmentsY_length * y_index, -0.5f + m_segmentsX_length * (x_index + 1), 1.0f }].insert({ -1.0f, 0.0f, 0.0f });
				normals_by_vertex[{ -0.5f, -0.5f + m_segmentsY_length * (y_index + 1), -0.5f + m_segmentsX_length * (x_index + 1), 1.0f }].insert({ -1.0f, 0.0f, 0.0f });

				normals_by_vertex[{ -0.5f, -0.5f + m_segmentsY_length * y_index, -0.5f + m_segmentsX_length * x_index, 1.0f }].insert({ -1.0f, 0.0f, 0.0f });
				normals_by_vertex[{ -0.5f, -0.5f + m_segmentsY_length * (y_index + 1), -0.5f + m_segmentsX_length * (x_index + 1), 1.0f }].insert({ -1.0f, 0.0f, 0.0f });
				normals_by_vertex[{ -0.5f, -0.5f + m_segmentsY_length * (y_index + 1), -0.5f + m_segmentsX_length * x_index, 1.0f }].insert({ -1.0f, 0.0f, 0.0f });
			}
		}

		vector<vector<float>> right_points = {};

		for (int x_index = 0; x_index < m_segmentsX; x_index++) {
			for (int y_index = 0; y_index < m_segmentsY; y_index++) {
				right_points.push_back({ 0.5f, -0.5f + m_segmentsY_length * y_index, -0.5f + m_segmentsX_length * (x_index + 1), 1.0f });
				right_points.push_back({ 0.5f, -0.5f + m_segmentsY_length * y_index, -0.5f + m_segmentsX_length * x_index, 1.0f });
				right_points.push_back({ 0.5f, -0.5f + m_segmentsY_length * (y_index + 1), -0.5f + m_segmentsX_length * x_index, 1.0f });

				right_points.push_back({ 0.5f, -0.5f + m_segmentsY_length * y_index, -0.5f + m_segmentsX_length * (x_index + 1), 1.0f });
				right_points.push_back({ 0.5f, -0.5f + m_segmentsY_length * (y_index + 1), -0.5f + m_segmentsX_length * x_index, 1.0f });
				right_points.push_back({ 0.5f, -0.5f + m_segmentsY_length * (y_index + 1), -0.5f + m_segmentsX_length * (x_index + 1), 1.0f });

				normals_by_vertex[{ 0.5f, -0.5f + m_segmentsY_length * y_index, -0.5f + m_segmentsX_length * (x_index + 1), 1.0f }].insert({ 1.0f, 0.0f, 0.0f });
				normals_by_vertex[{ 0.5f, -0.5f + m_segmentsY_length * y_index, -0.5f + m_segmentsX_length * x_index, 1.0f }].insert({ 1.0f, 0.0f, 0.0f });
				normals_by_vertex[{ 0.5f, -0.5f + m_segmentsY_length * (y_index + 1), -0.5f + m_segmentsX_length * x_index, 1.0f }].insert({ 1.0f, 0.0f, 0.0f });

				normals_by_vertex[{ 0.5f, -0.5f + m_segmentsY_length * y_index, -0.5f + m_segmentsX_length * (x_index + 1), 1.0f }].insert({ 1.0f, 0.0f, 0.0f });
				normals_by_vertex[{ 0.5f, -0.5f + m_segmentsY_length * (y_index + 1), -0.5f + m_segmentsX_length * x_index, 1.0f }].insert({ 1.0f, 0.0f, 0.0f });
				normals_by_vertex[{ 0.5f, -0.5f + m_segmentsY_length * (y_index + 1), -0.5f + m_segmentsX_length * (x_index + 1), 1.0f }].insert({ 1.0f, 0.0f, 0.0f });
			}
		}

		glBegin(GL_TRIANGLES);
		glNormal3f(0.0f, 0.0f, 1.0f);
		for (int index = 0; index < front_points.size(); index++) {
			glVertex3f(front_points[index][0], front_points[index][1], front_points[index][2]);
		}

		glNormal3f(0.0f, 0.0f, -1.0f);
		for (int index = 0; index < back_points.size(); index++) {
			glVertex3f(back_points[index][0], back_points[index][1], back_points[index][2]);
		}

		glNormal3f(0.0f, 1.0f, 0.0f);
		for (int index = 0; index < top_points.size(); index++) {
			glVertex3f(top_points[index][0], top_points[index][1], top_points[index][2]);
		}

		glNormal3f(0.0f, -1.0f, 0.0f);
		for (int index = 0; index < bottom_points.size(); index++) {
			glVertex3f(bottom_points[index][0], bottom_points[index][1], bottom_points[index][2]);
		}

		glNormal3f(-1.0f, 0.0f, 0.0f);
		for (int index = 0; index < left_points.size(); index++) {
			glVertex3f(left_points[index][0], left_points[index][1], left_points[index][2]);
		}

		glNormal3f(1.0f, 0.0f, 0.0f);
		for (int index = 0; index < right_points.size(); index++) {
			glVertex3f(right_points[index][0], right_points[index][1], right_points[index][2]);
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
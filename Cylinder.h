#ifndef CYLINDER_H
#define CYLINDER_H

#include "Shape.h"

class Cylinder : public Shape {
public:
	Cylinder() {};
	~Cylinder() {};
	std::vector<std::vector<float>> my_points;
	std::vector<std::vector<float>> my_normals;
	std::vector<std::vector<float>> my_normals_vector;

	OBJ_TYPE getType() {
		return SHAPE_CYLINDER;
	}


	void draw() {
		this->my_points.clear();
		this->my_normals.clear();
		this->my_normals_vector.clear();
		int segment_x = this->m_segmentsX;
		int segment_y = this->m_segmentsY;

		float x1 = 0;
		float y1 = 1;
		float z1 = 0;
		float w1 = 0;
		// if(my_points.empty()){
		vector<vector<float>> my_points_t = up_and_down(segment_x, segment_y);
		for (int p = 0; p < segment_y; p++) {
			my_points_t.push_back(my_points_t[p]);
			my_normals_vector.push_back(my_normals_vector[p]);
		}
		my_points = my_points_t;
		int index = 0;

		glBegin(GL_TRIANGLES);
		index = 0;
		for (int i = 0; i < segment_x; i++) {
			for (int p = 0; p < segment_y - 1; p++) {
				//render each square (made up of two triangles)
				int pnt_one_idx = i * segment_y + p;
				int pnt_two_idx = i * segment_y + p + 1;
				int pnt_three_idx = i * segment_y + p + segment_y;
				int pnt_four_idx = i * segment_y + p + 1 + segment_y;

				std::vector<float> val_1_n = my_normals_vector[pnt_one_idx];
				std::vector<float> val_1 = my_points[pnt_one_idx];
				std::vector<float> val_2_n = my_normals_vector[pnt_two_idx];
				std::vector<float> val_2 = my_points[pnt_two_idx];
				std::vector<float> val_3_n = my_normals_vector[pnt_three_idx];
				std::vector<float> val_3 = my_points[pnt_three_idx];
				std::vector<float> val_4_n = my_normals_vector[pnt_four_idx];
				std::vector<float> val_4 = my_points[pnt_four_idx];

				normalizeNormal(val_1_n[0], val_1_n[1], val_1_n[2]);
				glVertex3f(val_1[0], val_1[1], val_1[2]);
				normalizeNormal(val_3_n[0], val_3_n[1], val_3_n[2]);
				glVertex3f(val_3[0], val_3[1], val_3[2]);
				normalizeNormal(val_4_n[0], val_4_n[1], val_4_n[2]);
				glVertex3f(val_4[0], val_4[1], val_4[2]);

				normalizeNormal(val_1_n[0], val_1_n[1], val_1_n[2]);
				glVertex3f(val_1[0], val_1[1], val_1[2]);
				normalizeNormal(val_2_n[0], val_2_n[1], val_2_n[2]);
				glVertex3f(val_2[0], val_2[1], val_2[2]);
				normalizeNormal(val_4_n[0], val_4_n[1], val_4_n[2]);
				glVertex3f(val_4[0], val_4[1], val_4[2]);

			}
		}
		//make the bottom
		vector<float> top_cyl = { 0, .5, 0 };
		vector<float> bottom_cyl = { 0, -.5, 0 };
		for (int i = 0; i < segment_x; i++) {
			int pnt_one_idx = segment_y * i;
			int pnt_two_idx = segment_y * i + segment_y;
			vector<float> pnt_one = my_points[pnt_one_idx];
			vector<float> pnt_two = my_points[pnt_two_idx];
			//setNormal(pnt_one[0], pnt_one[1], pnt_one[2], pnt_two[0], pnt_two[1], pnt_two[2], bottom_cyl[0], bottom_cyl[1], bottom_cyl[2]);
			glVertex3f(pnt_one[0], pnt_one[1], pnt_one[2]);
			glVertex3f(pnt_two[0], pnt_two[1], pnt_two[2]);
			glVertex3f(bottom_cyl[0], bottom_cyl[1], bottom_cyl[2]);
		}


		for (int i = 0; i < segment_x; i++) {
			int pnt_one_idx = (segment_y * i + segment_y) - 1;
			int pnt_two_idx = (segment_y * i + segment_y + segment_y) - 1;
			vector<float> pnt_one = my_points[pnt_one_idx];
			vector<float> pnt_two = my_points[pnt_two_idx];
			//setNormal(pnt_one[0], pnt_one[1], pnt_one[2], top_cyl[0], top_cyl[1], top_cyl[2], pnt_two[0], pnt_two[1], pnt_two[2]);
			glVertex3f(pnt_one[0], pnt_one[1], pnt_one[2]);
			glVertex3f(top_cyl[0], top_cyl[1], top_cyl[2]);
			glVertex3f(pnt_two[0], pnt_two[1], pnt_two[2]);
		}
		glEnd();
	};

	void drawNormal() {
		glBegin(GL_LINES);
		for (int i = 0; i < this->my_normals.size(); i += 2) {
			vector<float> f_val = this->my_normals[i];
			vector<float> s_val = this->my_normals[i + 1];
			glVertex3f(f_val[0], f_val[1], f_val[2]);
			glVertex3f(s_val[0], s_val[1], s_val[2]);
		}
		glEnd();
	};


private:

	std::vector<std::vector<float>> up_and_down(int seg_x, int seg_y) {
		float x, y, z, r, z_step, normal_r;
		vector<vector<float>> my_vert_vals;
		normal_r = .6;
		float normal_x, normal_y;
		for (int j = 0; j < seg_x; j++) {
			z_step = 1.0 / (seg_y - 1);
			// cout << "z step is " << z_step << endl;
			// cout << "segments y is " << seg_y << endl;
			for (int i = 0; i < seg_y; i++) {
				r = .5;
				x = r * cos(2 * PI * j / seg_x);
				y = r * sin(2 * PI * j / seg_x);
				normal_x = normal_r * cos(2 * PI * j / seg_x);
				normal_y = normal_r * sin(2 * PI * j / seg_x);
				z = -.5 + z_step * i;
				// cout << "z is " << z << endl;

				// the z axis for us is how far up and down the point is (y in our case)
				my_vert_vals.push_back(vector<float>{x, z, y, (float)1.0});
				this->my_normals.push_back(vector<float>{x, z, y, (float)1.0});
				this->my_normals.push_back(vector<float>{normal_x, z, normal_y, (float)1.0});
				this->my_normals_vector.push_back(vector<float>{normal_x, z, normal_y, (float)1.0});

			}
		}
		return my_vert_vals;
	}


};

#endif
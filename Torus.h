#ifndef TORUS_H
#define TORUS_H

#include "Shape.h"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

class Torus : public Shape {
public:
	Torus() {};
	~Torus() {};

	OBJ_TYPE getType() {
		return SHAPE_SPECIAL1;
	}

	void draw() {
		glPushMatrix();
		glRotatef(90, 1, 0, 0);
		drawBody();
		glPopMatrix();
	};

	void drawNormal() {
		glPushMatrix();
		glRotatef(90, 1, 0, 0);
		drawBodyNormal();
		glPopMatrix();
	};

	float Intersect(glm::vec3 eyePointP, glm::vec3 rayV, glm::mat4 transformMatrix) {
		float t = MIN_ISECT_DISTANCE;
		return t;
	}

	glm::vec3 findIsectNormal(glm::vec3 eyePoint, glm::vec3 ray, float dist) {
		glm::vec3 normalVector(0.0f, 0.0f, 0.0f);
		return normalVector;
	}

private:

	void drawBody() {
		float radius = 0.15;
		float radiusRing = 0.5;
		float angle_delta = 2 * PI / (float)m_segmentsX;
		float angleRing_delta = 2 * PI / (float)m_segmentsY;

		float angleRing = 0;

		glBegin(GL_TRIANGLES);
		for (int i = 0; i < m_segmentsY; i++) {
			float angle = 0;
			for (int j = 0; j < m_segmentsX; j++) {
				glm::vec4 origin(0, 0, 0, 1.0f);
				glm::vec4 p (radius, 0, 0, 1.0f);
				glm::mat4 rotZmat = glm::rotate(glm::mat4(1.0f), angle, glm::vec3(0.0f, 0.0f, 1.0f));
				glm::vec4 pt = rotZmat * p;
				glm::mat4 rotZmat2 = glm::rotate(glm::mat4(1.0f), angle+angle_delta, glm::vec3(0.0f, 0.0f, 1.0f));
				glm::vec4 newPt = rotZmat2 * p;

				glm::mat4 translate = glm::translate(glm::mat4(1.0f), glm::vec3(radiusRing - radius, 0, 0));
				glm::mat4 rotY = glm::rotate(glm::mat4(1.0f), angleRing, glm::vec3(0, 1, 0));
				glm::mat4 transform1 = rotY * translate;

				glm::vec4 origin_ring = transform1 * origin;
				glm::vec4 pt_ring =  transform1 * pt;
				glm::vec4 newPt_ring = transform1 * newPt;

				glm::mat4 rotY2 = glm::rotate(glm::mat4(1.0f), angleRing + angleRing_delta, glm::vec3(0, 1, 0));
				glm::mat4 transform2 = rotY2 * translate;

				glm::vec4 origin_ring_next = transform2 * origin;
				glm::vec4 pt_ring_next =  transform2 * pt;
				glm::vec4 newPt_ring_next = transform2  * newPt;

				normalizeNormal(newPt_ring_next - origin_ring_next);
				glVertex3d(newPt_ring_next.x, newPt_ring_next.y, newPt_ring_next.z);
				normalizeNormal(newPt_ring - origin_ring);
				glVertex3d(newPt_ring.x, newPt_ring.y, newPt_ring.z);
				normalizeNormal(pt_ring - origin_ring);
				glVertex3d(pt_ring.x, pt_ring.y, pt_ring.z);

				normalizeNormal(pt_ring_next - origin_ring_next);
				glVertex3d(pt_ring_next.x, pt_ring_next.y, pt_ring_next.z);
				normalizeNormal(newPt_ring_next - origin_ring_next);
				glVertex3d(newPt_ring_next.x, newPt_ring_next.y, newPt_ring_next.z);
				normalizeNormal(pt_ring - origin_ring);
				glVertex3d(pt_ring.x, pt_ring.y, pt_ring.z);

				angle = angle + angle_delta;
			}
			angleRing = angleRing + angleRing_delta;
		}
		glEnd();
	};

	void drawBodyNormal() {

		float radius = 0.15;
		float radiusRing = 0.5;
		float angle_delta = 2 * PI / (float)m_segmentsX;
		float angleRing_delta = 2 * PI / (float)m_segmentsY;

		glBegin(GL_LINES);
		float angleRing = 0;
		for (int i = 0; i < m_segmentsY; i++) {
			float angle = 0;
			for (int j = 0; j < m_segmentsX; j++) {
				glm::vec4 origin(0, 0, 0, 1.0f);
				glm::vec4 p(radius, 0, 0, 1.0f);
				glm::mat4 rotZmat = glm::rotate(glm::mat4(1.0f), angle, glm::vec3(0.0f, 0.0f, 1.0f));
				glm::vec4 pt = rotZmat * p;
				glm::mat4 rotZmat2 = glm::rotate(glm::mat4(1.0f), angle + angle_delta, glm::vec3(0.0f, 0.0f, 1.0f));
				glm::vec4 newPt = rotZmat2 * p;

				glm::mat4 translate = glm::translate(glm::mat4(1.0f), glm::vec3(radiusRing - radius, 0, 0));
				glm::mat4 rotY = glm::rotate(glm::mat4(1.0f), angleRing, glm::vec3(0, 1, 0));
				glm::mat4 transform1 = rotY * translate;

				glm::vec4 origin_ring = transform1 * origin;
				glm::vec4 pt_ring = transform1 * pt;
				glm::vec4 newPt_ring = transform1 * newPt;

				glm::mat4 rotY2 = glm::rotate(glm::mat4(1.0f), angleRing + angleRing_delta, glm::vec3(0, 1, 0));
				glm::mat4 transform2 = rotY2 * translate;

				glm::vec4 origin_ring_next = transform2 * origin;
				glm::vec4 pt_ring_next = transform2 * pt;
				glm::vec4 newPt_ring_next = transform2 * newPt;

				glm::vec4 v = pt_ring - origin_ring;

				v = normalize(v);
				glVertex3d(pt_ring.x, pt_ring.y, pt_ring.z); glVertex4fv(glm::value_ptr(pt_ring+v*0.1f));
				v = newPt_ring - origin_ring;
				v = normalize(v);
				glVertex4fv(glm::value_ptr(newPt_ring)); glVertex4fv(glm::value_ptr(newPt_ring+v*0.1f)); 
				v = newPt_ring_next - origin_ring_next;
				v = normalize(v);
				glVertex4fv(glm::value_ptr(newPt_ring_next)); glVertex4fv(glm::value_ptr(newPt_ring_next+v*0.1f)); 

				v = (newPt_ring_next - origin_ring_next);
				v = normalize(v);
				glVertex4fv(glm::value_ptr(newPt_ring_next)); glVertex4fv(glm::value_ptr(newPt_ring_next + v * 0.1f));
				v = (pt_ring_next - origin_ring_next);
				v= normalize(v);
				glVertex4fv(glm::value_ptr(pt_ring_next)); glVertex4fv(glm::value_ptr(pt_ring_next+v*0.1f)); 
				v = (pt_ring - origin_ring);
				v = normalize(v);
				glVertex4fv(glm::value_ptr(pt_ring)); glVertex4fv(glm::value_ptr(pt_ring+v*0.1f)); 

				angle = angle + angle_delta;
			}
			angleRing = angleRing + angleRing_delta;
		}
		glEnd();
	};
};

#endif
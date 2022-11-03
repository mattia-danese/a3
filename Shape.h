#ifndef SHAPE_H
#define SHAPE_H

#include <FL/gl.h>
#include <FL/glu.h>
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <vector>

#define PI glm::pi<float>()

enum OBJ_TYPE {
	SHAPE_CUBE = 0,
	SHAPE_CYLINDER = 1,
	SHAPE_CONE = 2,
	SHAPE_SPHERE = 3,
	SHAPE_SPECIAL1 = 4,
	SHAPE_SPECIAL2 = 5,
	SHAPE_SPECIAL3 = 6,
	SHAPE_MESH = 7,
	SHAPE_SCENE = 8
};

class Shape {
public:
	static int m_segmentsX;
	static int m_segmentsY;

	Shape() {
		setSegments(10, 10);
	};
	~Shape() {};

	void setSegments(int segX, int segY) {
		m_segmentsX = segX;
		m_segmentsY = segY;
	}

	virtual OBJ_TYPE getType() = 0;
	virtual void draw() {};
	virtual void drawNormal() {};

protected:
	void normalizeNormal(float x, float y, float z) {
		normalizeNormal(glm::vec3(x, y, z));
	};

	void normalizeNormal(glm::vec3 v) {
		glm::vec3 tmpV = glm::normalize(v);
		glNormal3f(tmpV.x, tmpV.y, tmpV.z);
	};

    std::vector<float> getNormal(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3) {
        float v1x, v1y, v1z;
        float v2x, v2y, v2z;
        float cx, cy, cz;

        // find vector between x2 and x1
        v1x = x1 - x2;
        v1y = y1 - y2;
        v1z = z1 - z2;

        // find vector between x3 and x2
        v2x = x2 - x3;
        v2y = y2 - y3;
        v2z = z2 - z3;

        // cross product v1xv2
        cx = v1y * v2z - v1z * v2y;
        cy = v1z * v2x - v1x * v2z;
        cz = v1x * v2y - v1y * v2x;

        // normalize
        float length = sqrt(cx * cx + cy * cy + cz * cz);
        cx = cx / length;
        cy = cy / length;
        cz = cz / length;

        return { cx, cy, cz };
    };

    std::vector<float> getNormalizedVector(float cx, float cy, float cz) {
        // normalize
        float length = sqrt(cx * cx + cy * cy + cz * cz);
        cx = cx / length;
        cy = cy / length;
        cz = cz / length;

        return { cx, cy, cz };
    }

    void setNormal(float x1, float y1, float z1, float x2, float y2, float z2,
        float x3, float y3, float z3) {

        float v1x, v1y, v1z;
        float v2x, v2y, v2z;
        float cx, cy, cz;

        // find vector between x2 and x1
        v1x = x1 - x2;
        v1y = y1 - y2;
        v1z = z1 - z2;

        // find vector between x3 and x2
        v2x = x2 - x3;
        v2y = y2 - y3;
        v2z = z2 - z3;

        // cross product v1xv2

        cx = v1y * v2z - v1z * v2y;
        cy = v1z * v2x - v1x * v2z;
        cz = v1x * v2y - v1y * v2x;

        // normalize

        float length = sqrt(cx * cx + cy * cy + cz * cz);
        cx = cx / length;
        cy = cy / length;
        cz = cz / length;

        glNormal3f(cx, cy, cz);
    }
};

#endif
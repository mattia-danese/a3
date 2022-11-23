#pragma once

#ifndef MYGLCANVAS_H
#define MYGLCANVAS_H

#include <FL/gl.h>
#include <FL/glut.h>
#include <FL/glu.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <time.h>
#include <iostream>
#include "list"

#include "Shape.h"
#include "Cube.h"
#include "Cylinder.h"
#include "Cone.h"
#include "Sphere.h"
#include "ppm.h"

#include "Camera.h"
#include "scene/SceneParser.h"

class MyGLCanvas : public Fl_Gl_Window {
public:

	struct nearestObj {
				float t;
				ScenePrimitive* prim; 
				glm::mat4 m_inv;
	};

	glm::vec3 rotVec;
	glm::vec3 eyePosition;
	GLubyte* pixels = NULL;

	int isectOnly;
	int segmentsX, segmentsY;
	int depth;
	int depth_fresh;
	int samples_per_pixel;
	float scale;

	OBJ_TYPE objType;
	Cube* cube;
	Cylinder* cylinder;
	Cone* cone;
	Sphere* sphere;
	Shape* shape;

	Camera* camera;
	SceneParser* parser;

	// what I added
	list<ScenePrimitive*> primitives;
	vector<SceneTransformation*> scenetransformations;
	vector<pair<ScenePrimitive*, vector<SceneTransformation*>>> my_scene_vals;
	vector<vector<int>> isect_pixels;
	glm::vec3 lookAt;
	SceneCameraData cameraData;
	map<string, ppm*> p_map;
	glm::vec3 ist_min;
	ScenePrimitive* prim_m;
	ScenePrimitive* original_prim;
	glm::vec3 cube_n;

	MyGLCanvas(int x, int y, int w, int h, const char *l = 0);
	~MyGLCanvas();
	void renderShape(OBJ_TYPE type);
	void setSegments();
	void loadSceneFile(const char* filenamePath);
	void renderScene();
	SceneColor computeColor(ScenePrimitive* p_m, glm::vec3 Nhat, glm::vec3 pos);
	glm::vec3 computeNormal(glm::vec3 inst, OBJ_TYPE shape);
	void traverse1(SceneNode* root, vector<pair<ScenePrimitive*, vector<SceneTransformation*>>>& my_scene_vals, vector<SceneTransformation*> curr_trans,map<string, ppm*>* p_m);
	float intersectCube (glm::vec3 eyePointP, glm::vec3 rayV, glm::mat4 transformMatrix, glm::vec3 cubepos);
	double intersectSphere(glm::vec3 eyePointP, glm::vec3 rayV, glm::mat4 transformMatrix, glm::vec3 spherepos);
	double intersectCone(glm::vec3 eyePointP, glm::vec3 rayV, glm::mat4 transformMatrix, glm::vec3 spherepos);
	double intersectCylinder (glm::vec3 eyePointP, glm::vec3 rayV, glm::mat4 transformMatrix, glm::vec3 spherepos);
	float quadraticForm(float A, float B, float C);
	void getst(int* s, int* t, int u, int v, int w, int h, int i, int j);
// 	double intersectCube (glm::vec3 eyePointP, glm::vec3 rayV, glm::mat4 transformMatrix, glm::vec3 spherepos);
	float intersectsq(glm::vec3 eye, glm::vec3 d, int i, float n);
	glm::vec3 camLookAt(int pixelX, int pixelY);
	SceneColor bound(SceneColor c);
	bool shadowCheck(glm::vec3 pos, glm::vec3 Lhati);
	SceneColor textureMap(SceneColor color, ScenePrimitive* prim);
	void convert_xyz_to_cube_uv(float x, float y, float z, int *index, float *u, float *v);
	SceneColor loopObjects(vector<pair<ScenePrimitive*, vector<SceneTransformation*>>> my_scene_vals, glm::vec3 eye_pnt, glm::vec3 ray, int* hit, bool shadow_check);
	float mPos(float x, float y);

private:

	void setpixel(GLubyte* buf, int x, int y, int r, int g, int b);

	glm::vec3 generateRay(int pixelX, int pixelY);
	glm::vec3 getEyePoint();
	glm::vec3 getIsectPointWorldCoord(glm::vec3 eye, glm::vec3 ray, float t);

	void draw();

	int handle(int);
	void resize(int x, int y, int w, int h);
	void updateCamera(int width, int height);

	int pixelWidth, pixelHeight;
};

#endif // !MYGLCANVAS_H

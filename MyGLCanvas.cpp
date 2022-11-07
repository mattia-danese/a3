#define NUM_OPENGL_LIGHTS 8

#include "MyGLCanvas.h"
#include <glm/gtx/string_cast.hpp>
#include <float.h>

int Shape::m_segmentsX;
int Shape::m_segmentsY;

MyGLCanvas::MyGLCanvas(int x, int y, int w, int h, const char *l) : Fl_Gl_Window(x, y, w, h, l) {
	mode(FL_RGB | FL_ALPHA | FL_DEPTH | FL_DOUBLE);
	
	rotVec = glm::vec3(0.0f, 0.0f, 0.0f);
	eyePosition = glm::vec3(2.0f, 2.0f, 2.0f);

	pixelWidth = w;
	pixelHeight = h;

	isectOnly = 1;
	segmentsX = segmentsY = 10;
	scale = 1.0f;
	parser = NULL;

	objType = SHAPE_CUBE;
	cube = new Cube();
	cylinder = new Cylinder();
	cone = new Cone();
	sphere = new Sphere();
	shape = cube;

	shape->setSegments(segmentsX, segmentsY);
	camera = new Camera();
	camera->orientLookVec(eyePosition, glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));
}

MyGLCanvas::~MyGLCanvas() {
	delete cube;
	delete cylinder;
	delete cone;
	delete sphere;
	if (camera != NULL) {
		delete camera;
	}
	if (parser != NULL) {
		delete parser;
	}
	if (pixels != NULL) {
		delete pixels;
	}
}

void MyGLCanvas::renderShape(OBJ_TYPE type) {
	objType = type;
	switch (type) {
	case SHAPE_CUBE:
		shape = cube;
		break;
	case SHAPE_CYLINDER:
		shape = cylinder;
		break;
	case SHAPE_CONE:
		shape = cone;
		break;
	case SHAPE_SPHERE:
		shape = sphere;
		break;
	case SHAPE_SPECIAL1:
	default:
		shape = cube;
	}

	shape->setSegments(segmentsX, segmentsY);
	shape->draw();
}

void MyGLCanvas::loadSceneFile(const char* filenamePath) {
	my_scene_vals.clear();
	if (parser != NULL) {
		delete parser;
	}
	parser = new SceneParser(filenamePath);

	bool success = parser->parse();
	std::cout << "success? " << success << endl;
	if (success == false) {
		delete parser;
		parser = NULL;
	}
	else {
		SceneCameraData cameraData;
		parser->getCameraData(cameraData);

		camera->reset();
		camera->setViewAngle(cameraData.heightAngle);
		if (cameraData.isDir == true) {
			camera->orientLookVec(cameraData.pos, cameraData.look, cameraData.up);
		}
		else {
			camera->orientLookAt(cameraData.pos, cameraData.lookAt, cameraData.up);
		}
	}
}


void MyGLCanvas::setSegments() {
	shape->setSegments(segmentsX, segmentsY);
}

glm::vec3 MyGLCanvas::getEyePoint() {
	return camera->getEyePoint();
	// glm::mat4 inv =  camera->getInverseModelViewMatrix() * camera->getInverseScaleMatrix();
	// return glm::vec3(inv * glm::vec4(camera->getEyePoint(), 1.0));
}

glm::vec3 MyGLCanvas::generateRay(int pixelX, int pixelY) {
	glm::vec3 lookAt = glm::vec3(-1.0f + 2.0f * ((float)pixelX / (float)camera->getScreenWidth()),
		-1.0f + 2.0f * ((float)pixelY / (float)camera->getScreenHeight()),
		-1.0f);

	//likely stems from here
	glm::mat4 inv =  camera->getInverseModelViewMatrix() * camera->getInverseScaleMatrix();
	glm::vec4 worldFPP = inv * glm::vec4(lookAt, 0.0f);
	// already in world coordinates
	glm::vec4 worldEye = glm::vec4(getEyePoint(), 1.0);

	glm::vec4 dHat = glm::normalize(worldFPP - worldEye);

	// std::cout << "X: " << dHat.x << "," << "Y: " << dHat.y << ',' << "Z: " << dHat.z << std::endl;

	return glm::vec3(dHat.x, dHat.y, dHat.z);
}

glm::vec3 MyGLCanvas::getIsectPointWorldCoord(glm::vec3 eye, glm::vec3 ray, float t) {
	glm::vec3 p = eye + (t * ray);
	return p;
}

// FROM PREVIOUS LAB
double MyGLCanvas::intersectSphere (glm::vec3 eyePointP, glm::vec3 rayV, glm::mat4 transformMatrix, glm::vec3 spherepos) {
	glm::mat4 transformInv = glm::inverse(transformMatrix);
	glm::vec3 d = glm::vec3(transformInv * glm::vec4(rayV,0));
	glm::vec3 eye = glm::vec3(transformInv * glm::vec4(eyePointP, 1));
	// std::cout << glm::to_string(transformInv) << std::endl;
	// glm::vec3 eyePointObject = glm::vec3(transformInv * glm::vec4(eyePointP, 1));
    // glm::vec3 rayVObject = glm::vec3(transformInv * glm::vec4(rayV, 0));

	double A = glm::dot(d, d);
	double B = 2 * glm::dot(eye, d);
	double C = glm::dot(eye, eye) - pow(0.5, 2);
	double discriminant = pow(B, 2) - (4.0 * A * C);

    // std::cout << "Ray: " << rayV.x << "," << rayV.y << "," << rayV.z << std::endl;
    // std::cout << "Eye: " << eyePointP.x << "," << eyePointP.y << "," << eyePointP.z << std::endl;

    // std::cout << "A: " << A << std::endl;
    // std::cout << "B: " << B << std::endl;
    // std::cout << "C: " << C << std::endl;
    // std::cout << "discriminant: " << discriminant << std::endl;

	if (discriminant < 0) {return -1;}
    else {
        double t1 = (-B + sqrt(discriminant)) / (2.0 * A);
        double t2 = (-B - sqrt(discriminant)) / (2.0 * A);

        if (t1 < 0) {return t2;}
        if (t2 < 0) {return t1;}
        return std::min(t1, t2);
    }
}

float MyGLCanvas::solveQuadratic(double A, double B, double C) {
    double det = B*B - 4.0*A*C;
    //std::cout << "det is  " << det << std::endl;
    if (det > 0) {
        double t0 = (-1.0*B + sqrt(det)) / (2.0*A);
        double t1 = (-1.0*B - sqrt(det)) / (2.0*A);
        // std::cout << "t0 is  " << t0 << std::endl;
		// std::cout << "t1 is  " << t1 << std::endl;
        return min(t0, t1);
    } else if (det == 0) {
        double t = (-B + sqrt(det)) / (2.0*A);

        return t;
    } else {
        return -1; 
    }   
}


bool inBounds(double p1, double p2) {
	if (-0.5 <= p1 and p1 <= 0.5) {
		if (-0.5 <= p2 and p2 <= 0.5)
			return true;
	}
	return false;
}


float MyGLCanvas::intersectCube (glm::vec3 eyePointP, glm::vec3 rayV, glm::mat4 transformMatrix, glm::vec3 spherepos) {
	glm::mat4 transformInv = glm::inverse(transformMatrix);
	glm::vec4 d = transformInv * glm::vec4(rayV,0);
	glm::vec4 eye = transformInv * glm::vec4(eyePointP, 1);

	float t1 = intersectFace(eye, d, 0, 0.5);
    float t2 = intersectFace(eye, d, 1, 0.5);
    float t3 = intersectFace(eye, d, 2, 0.5);
    float t4 = intersectFace(eye, d, 0, -0.5);
    float t5 = intersectFace(eye, d, 1, -0.5);
    float t6 = intersectFace(eye, d, 2, -0.5);
    
    return minPos(t1, (minPos(t2, (minPos(t3, (minPos(t4, (minPos(t5, t6)))))))));
}

float MyGLCanvas::intersectFace(glm::vec3 eye, glm::vec3 d, int i, float n) {
    float t = (n - eye[i]) / d[i];
    glm::vec3 intersect = eye + d * t;
    if ((intersect[(i + 1) % 3] < 0.5 && intersect[(i + 1) % 3] > -0.5) &&
        (intersect[(i + 2) % 3] < 0.5 && intersect[(i + 2) % 3] > -0.5)) {
        return t;
    } else {
        return -1;
    }
}

double MyGLCanvas::intersectCone (glm::vec3 eyePointP, glm::vec3 rayV, glm::mat4 transformMatrix, glm::vec3 spherepos) {
	glm::mat4 transformInv = glm::inverse(transformMatrix);
	glm::vec3 eyePointObject = transformInv * glm::vec4(eyePointP, 1);
    glm::vec3 rayVObject = transformInv * glm::vec4(rayV, 0);
	glm::vec3 d = rayVObject;
	glm::vec3 eye = eyePointObject;
    double A = d[0]*d[0] + d[2]*d[2] - (.25*d[1]*d[1]);
    double B = 2*eye[0]*d[0] + 2*eye[2]*d[2] - .5*eye[1]*d[1] + .25*d[1];
    double C = eye[0]*eye[0] + eye[2]*eye[2] - .25*eye[1]*eye[1] + .25*eye[1] - .25*.25;

    float t_sides = solveQuadratic(A, B, C);
    //     std::cout << "Ray: " << rayV.x << "," << rayV.y << "," << rayV.z << std::endl;
    // std::cout << "Eye: " << eyePointP.x << "," << eyePointP.y << "," << eyePointP.z << std::endl;

    // std::cout << "A: " << A << std::endl;
    // std::cout << "B: " << B << std::endl;
    // std::cout << "C: " << C << std::endl;
    glm::vec3 intersect = eye +  d *t_sides;
	if (!(intersect[1] > -0.5 && intersect[1] < 0.5)) {
        t_sides = -1;
    }

    float t_cap = (-0.5 - eye[1]) / d[1];
    intersect = eye + d * t_cap;
    if (!(intersect[0]*intersect[0] + intersect[2]*intersect[2] <= 0.25)) {
        t_cap = -1;
    }
    
    return minPos(t_sides, t_cap);
}

float MyGLCanvas::minPos(float x, float y) {
    if (x < 0) {
        if (y < 0)
            return -1;
        else 
            return y;
    } else {
        if (y < 0)
            return x;
        else 
            return min(x,y);
    }
}


//
//
// NEEDS IMPLEMENTING
//
//
double MyGLCanvas::intersectCylinder (glm::vec3 eyePointP, glm::vec3 rayV, glm::mat4 transformMatrix, glm::vec3 spherepos) {
	glm::mat4 transformInv = glm::inverse(transformMatrix);
	glm::vec3 eyePointObject = transformInv * glm::vec4(eyePointP, 1);
    glm::vec3 rayVObject = transformInv * glm::vec4(rayV, 0);
	glm::vec3 d = rayVObject;
	glm::vec3 eye = eyePointObject;
    double A = d[0]*d[0] + d[2]*d[2];
    double B = 2*eye[0]*d[0] + 2*eye[2]*d[2];
    double C = eye[0]*eye[0] + eye[2]*eye[2] - 0.25;


    // std::cout << "Ray: " << rayV.x << "," << rayV.y << "," << rayV.z << std::endl;
    // std::cout << "Eye: " << eyePointP.x << "," << eyePointP.y << "," << eyePointP.z << std::endl;

    // std::cout << "A: " << A << std::endl;
    // std::cout << "B: " << B << std::endl;
    // std::cout << "C: " << C << std::endl;

    float t_body = solveQuadratic(A, B, C);

	// return t_body;

	    // std::cout << "discriminant: " << discriminant << std::endl;
    glm::vec3 intersect = eye + d * t_body;
    if (!(intersect[1] > -0.5 && intersect[1] < 0.5)) {
        t_body = -1;
    }

    float t_cap1 = (0.5 - eye[1]) / d[1];
    intersect = eye + d * t_cap1;
    if (!(intersect[0]*intersect[0] + intersect[2]*intersect[2] <= 0.25)) {
        t_cap1 = -1;
    }

    float t_cap2 = (-0.5 - eye[1]) / d[1];
    intersect = eye + d * t_cap2;
    if (!(intersect[0]*intersect[0] + intersect[2]*intersect[2] <= 0.25)) {
        t_cap2 = -1;
    }
	// float t = t_body;
	// if(t_cap1 > 0 && t_cap2 > 0){
	// 	return min(t_cap2, t_cap1);
	// }else if(t_cap1 < 0 && t_cap2 >0){
	// 	t = t_cap2;
	// 	return t;
	// }else if(t_cap2 < 0 && t_cap1 >0){
	// 	t = t_cap1;
	// 	return t;
	// }
    return minPos(t_body, minPos(t_cap2, t_cap1));
}



void MyGLCanvas::draw() {
	if (!valid()) {  //this is called when the GL canvas is set up for the first time or when it is resized...
		printf("establishing GL context\n");

		glViewport(0, 0, w(), h());
		updateCamera(w(), h());

		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	}
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if (parser == NULL) {
		return;
	}

	if (pixels == NULL) {
		return;
	}

	//this just draws the "pixels" to the screen
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glDrawPixels(pixelWidth, pixelHeight, GL_RGB, GL_UNSIGNED_BYTE, pixels);
}

int MyGLCanvas::handle(int e) {
	//printf("Event was %s (%d)\n", fl_eventnames[e], e);
	switch (e) {
	case FL_KEYUP:
		printf("keyboard event: key pressed: %c\n", Fl::event_key());
		break;
	case FL_MOUSEWHEEL:
		break;
	}

	return Fl_Gl_Window::handle(e);
}

void MyGLCanvas::resize(int x, int y, int w, int h) {
	Fl_Gl_Window::resize(x, y, w, h);
	if (camera != NULL) {
		camera->setScreenSize(w, h);
	}
	puts("resize called");
}


void MyGLCanvas::updateCamera(int width, int height) {
	float xy_aspect;
	xy_aspect = (float)width / (float)height;
	camera->setScreenSize(width, height);
}

//Given the pixel (x, y) position, set its color to (r, g, b)
void MyGLCanvas::setpixel(GLubyte* buf, int x, int y, int r, int g, int b) {
	pixelWidth = camera->getScreenWidth();
	buf[(y*pixelWidth + x) * 3 + 0] = (GLubyte)r;
	buf[(y*pixelWidth + x) * 3 + 1] = (GLubyte)g;
	buf[(y*pixelWidth + x) * 3 + 2] = (GLubyte)b;
}

SceneColor MyGLCanvas::computeSceneColor(SceneMaterial material, glm::vec3 Nhat, glm::vec3 pos) {

	SceneGlobalData data;
	parser->getGlobalData(data);

	float ka = data.ka;

	float kd = data.kd;

	SceneColor Oa = material.cAmbient;
	SceneColor Od = material.cDiffuse;

	SceneColor color;
	color.r = 0;
	color.g = 0;
	color.b = 0;
	color.a = 255;

    int nLights = parser->getNumLights();
    for (int m = 0; m < nLights; m++) {
        SceneLightData lData;
        parser->getLightData(m, lData);
        SceneColor lm = lData.color;
        glm::vec3 Lhatm = lData.pos - pos;
        glm::normalize(Lhatm);
//         std::cout << "nhat is  " << Nhat[0] << Nhat[1] << Nhat[2] << std::endl;
        //std::cout << lm.g << std::endl;
       	if (Od.r > 1.0) Od.r = 1.0;
       	if (Od.g > 1.0) Od.g = 1.0;
       	if (Od.b > 1.0) Od.b = 1.0;
       	if (lm.r > 1.0) lm.r = 1.0;
       	if (lm.g > 1.0) lm.g = 1.0;
       	if (lm.b > 1.0) lm.b = 1.0;

        if (dot(Nhat, Lhatm) > 0) {
            color.r += (kd * Od.r * lm.r * dot(Nhat, Lhatm));  
            color.g += (kd * Od.g * lm.g * dot(Nhat, Lhatm));  
            color.b += (kd * Od.b * lm.b * dot(Nhat, Lhatm));  
        }
    }
    if (Oa.r > 1.0) Oa.r = 1.0;
   	if (Oa.g > 1.0) Oa.g = 1.0;
   	if (Oa.b > 1.0) Oa.b = 1.0;
    color.r += ka * (Oa.r);
    color.g += ka * (Oa.g);
    color.b += ka * (Oa.b);

    color.r *= 255;
    color.g *= 255;
    color.b *= 255;

	return color;
}

glm::vec3 MyGLCanvas::computeNormal(glm::vec3 intersection, OBJ_TYPE shape) {
	    switch(shape) {
        case SHAPE_CUBE: {
            if (intersection[0] > .4999)
                return glm::vec3(1, 0, 0);
            if (intersection[0] < -.4999)
                return glm::vec3(-1, 0, 0);
            if (intersection[1] > .4999)
                return glm::vec3(0, 1, 0);
            if (intersection[1] < -.4999)
                return glm::vec3(0, -1, 0);
            if (intersection[2] > .4999)
                return glm::vec3(0, 0, 1);
            if (intersection[2] < -.4999)
                return glm::vec3(0, 0, -1);
            }
            break;
        case SHAPE_CYLINDER:
           if (IN_RANGE(intersection[1], 0.5))
               return glm::vec3(0, 1, 0);
           if (IN_RANGE(intersection[1], -0.5))
               return glm::vec3(0, -1, 0);
           return glm::vec3(intersection[0], 0, intersection[2]);
           break; 
        case SHAPE_CONE: {
           if (IN_RANGE(intersection[1], -0.5))
               return glm::vec3(0, -1, 0);
           glm::vec3 v1 = glm::vec3(intersection[0], 0, intersection[2]);
           glm::normalize(v1);
           glm::vec3 v2 = glm::vec3(0, .5, 0);
           return v1 + v2; }
           break; 
        case SHAPE_SPHERE:
           return glm::vec3(intersection[0], intersection[1], intersection[2]);
           break;
    }
}

void MyGLCanvas::traverse1(SceneNode* root, vector<pair<ScenePrimitive*, vector<SceneTransformation*>>>& my_scene_vals, vector<SceneTransformation*> curr_trans)
{
	if (root == nullptr)
		return;

	for (SceneNode* node : root->children) { // storing transformations in non-primitive nodes 
		for (SceneTransformation* my_trans : node->transformations) {
			curr_trans.push_back(my_trans);
		}

		traverse1(node, my_scene_vals, curr_trans);

		// hit leaf node
		// for (SceneTransformation* my_trans : node->transformations) { //storing transformations in primitive node
		// 	curr_trans.push_back(my_trans);
		// }

		for (ScenePrimitive* my_prim : node->primitives) { // for all primitives in leaf node, make pair of prim and previously collected trasnformations
			list<SceneTransformation*> curr_scene_trans = {};
			my_scene_vals.push_back(pair<ScenePrimitive*, vector<SceneTransformation*>>(my_prim, curr_trans));
		}
		for (SceneTransformation* my_trans : node->transformations) { // func will recurse back up, delete transformations
			curr_trans.pop_back();
		}
	}
}


void MyGLCanvas::renderScene() {
	std::cout << "render button clicked!" << endl;

	if (parser == NULL) {
		std::cout << "no scene loaded yet" << endl;
		return;
	}

	pixelWidth = w();
	pixelHeight = h();


	updateCamera(pixelWidth, pixelHeight);

	if (pixels != NULL) {
		delete pixels;
	}
	pixels = new GLubyte[pixelWidth  * pixelHeight * 3];
	memset(pixels, 0, pixelWidth  * pixelHeight * 3);
		glm::vec3  eye_pnt = getEyePoint();
		if (my_scene_vals.empty()) {
			std::cout << "calling traverse!" << endl;
			SceneNode* root = parser->getRootNode();
			traverse1(root, my_scene_vals, scenetransformations);
		}

		for (int i = 0; i < pixelWidth; i++) {
			for (int j = 0; j < pixelHeight; j++) {
                glm::vec3 intersection_obj = glm::vec4(0);
				glm::vec4 intersection = glm::vec4(0);
				float t_min = FLT_MAX;
				SceneColor color;
				color.r = 0, color.g = 0, color.b = 0;
				//TODO: this is where your ray casting will happen!
				glm::vec3 ray = generateRay(i, j);
				for (pair<ScenePrimitive*, vector<SceneTransformation*>> my_p : my_scene_vals) {
					ScenePrimitive* prim = my_p.first;
					float t = -1;
					glm::mat4 m = glm::mat4(1.0);
					// loop through transformations on each prim creating the matrix
					for (SceneTransformation* my_trans : my_p.second) {
						switch (my_trans->type) {
							case(TRANSFORMATION_TRANSLATE):
								//std::cout << my_trans->translate.x << " " << my_trans->translate.y << " " << my_trans->translate.z << std::endl;
								m = m * glm::translate(glm::mat4(1.0), my_trans->translate);
								break;
							case(TRANSFORMATION_SCALE):
								//std::cout << my_trans->scale.x << " " << my_trans->scale.y << " " << my_trans->scale.z << std::endl;
								m = m * glm::scale(glm::mat4(1.0), my_trans->scale);
								break;
							case(TRANSFORMATION_ROTATE):
								m = m * glm::rotate(glm::mat4(1.0), my_trans->angle, my_trans->rotate);
								break;
							case(TRANSFORMATION_MATRIX):
								m = m * my_trans->matrix;
								break;
						}
					}
					glm::vec3 center = glm::vec3(0.0f);
					switch (prim->type) {
					case SHAPE_CUBE:
					//if (i == pixelWidth / 2 && j == pixelHeight / 2) {
						t = intersectCube(eye_pnt, ray, m, center);
					// 	std::cout << "eye point is " << eye_pnt.x << " " << eye_pnt.y << " " << eye_pnt.z << " " << std::endl;
					// 	std::cout << "ray is " << ray.x << " " << ray.y << " " << ray.z << " " << std::endl;
					// 	std::cout << "t is " << t << std::endl;
					// }
						break;
					case SHAPE_CYLINDER:
					//if (i == pixelWidth / 2 && j == pixelHeight / 2) {
						t = intersectCylinder(eye_pnt, ray, m, center);
					// 		std::cout << "eye point is " << eye_pnt.x << " " << eye_pnt.y << " " << eye_pnt.z << " " << std::endl;
					// 	std::cout << "ray is " << ray.x << " " << ray.y << " " << ray.z << " " << std::endl;
					// 	std::cout << "t is " << t << std::endl;
					// }
						break;
					case SHAPE_CONE:
						// cout << "here " << endl;
						// std::cout << "eye point is " << eye_pnt.x << " " << eye_pnt.y << " " << eye_pnt.z << " " << std::endl;
						// std::cout << "ray is " << ray.x << " " << ray.y << " " << ray.z << " " << std::endl;
						// if (i == pixelWidth / 2 && j == pixelHeight / 2) {
						t = intersectCone(eye_pnt, ray, m, center);
						// }
						break;
					case SHAPE_SPHERE:
							t = intersectSphere(eye_pnt, ray, m, center);
						break;
					}

					if (i == pixelWidth / 2 && j == pixelHeight / 2) {
						std::cout << "t is " << t << std::endl;
					}
				
					if (t > 0 && t < t_min) {
						// std::cout << "IGGGIGI" << endl;
						t_min = t;
						intersection_obj = getIsectPointWorldCoord(glm::vec3(glm::vec4(eye_pnt,1.0) * glm::inverse(m)), glm::vec3(glm::vec4(ray,0) * glm::inverse(m)), t);
						glm::vec3 normal = computeNormal(intersection_obj, prim->type); 
						glm::vec3 intersection = glm::vec3(transpose(m) * glm::vec4(intersection_obj, 1));
						normal = glm::vec3(glm::transpose(glm::inverse(m)) * glm::vec4(normal,0));
						glm::normalize(normal);
						color = computeSceneColor(prim->material, normal, intersection);
					}

				}

				if (isectOnly == 1) {
					setpixel(pixels, i, j, color.r, color.g, color.b);
				}
				else {
					setpixel(pixels, i, j, 255, 255, 255);
				}
			}
		}
	std::cout << "render complete" << endl;
	redraw();
}

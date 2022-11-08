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

	double A = glm::dot(d, d);
	double B = 2 * glm::dot(eye, d);
	double C = glm::dot(eye, eye) - pow(0.5, 2);
	double discriminant = pow(B, 2) - (4.0 * A * C);


	if (discriminant < 0) {return -1;}
    else {
        double t1 = (-B + sqrt(discriminant)) / (2.0 * A);
        double t2 = (-B - sqrt(discriminant)) / (2.0 * A);

        if (t1 < 0) {return t2;}
        if (t2 < 0) {return t1;}
        return std::min(t1, t2);
    }
}

float MyGLCanvas::quadraticForm(double A, double B, double C) {
    double det = B*B - 4.0*A*C;
    if (det > 0) {
        double t1 = (-1.0*B + sqrt(det)) / (2.0*A);
        double t2 = (-1.0*B - sqrt(det)) / (2.0*A);
		if (t1 < 0) {return t2;}
        if (t2 < 0) {return t1;}
        return std::min(t1, t2);
    } else {
        return -1; 
    }   
}

float MyGLCanvas::intersectCube (glm::vec3 eyePointP, glm::vec3 rayV, glm::mat4 transformMatrix, glm::vec3 spherepos) {
	glm::mat4 transformInv = glm::inverse(transformMatrix);
	glm::vec4 d_v = transformInv * glm::vec4(rayV,0);
	glm::vec4 eye = transformInv * glm::vec4(eyePointP, 1);

	float t4 = intersectsq(eye, d_v, 0, -0.5);
    float t5 = intersectsq(eye, d_v, 1, -0.5);
    float t6 = intersectsq(eye, d_v, 2, -0.5);

	float t1 = intersectsq(eye, d_v, 0, 0.5);
    float t2 = intersectsq(eye, d_v, 1, 0.5);
    float t3 = intersectsq(eye, d_v, 2, 0.5);
    
	//which face did we intersect
    return mPos(t1, (mPos(t2, (mPos(t3, (mPos(t4, (mPos(t5, t6)))))))));
}

float MyGLCanvas::intersectsq(glm::vec3 eye, glm::vec3 d, int i, float n) {
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
	glm::vec3 d_v = rayVObject;
	glm::vec3 eye_o = eyePointObject;
    double A = d_v[0]*d_v[0] + d_v[2]*d_v[2] - (.25*d_v[1]*d_v[1]);
    double B = 2*eye_o[0]*d_v[0] + 2*eye_o[2]*d_v[2] - .5*eye_o[1]*d_v[1] + .25*d_v[1];
    double C = eye_o[0]*eye_o[0] + eye_o[2]*eye_o[2] - .25*eye_o[1]*eye_o[1] + .25*eye_o[1] - .25*.25;

    float t_s = quadraticForm(A, B, C);
	float t_c = (-0.5 - eye_o[1]) / d_v[1];
    glm::vec3 intersect = eye_o +  d_v *t_s;
	if (!(intersect[1] > -0.5 && intersect[1] < 0.5)) {
        t_s = -1;
    }
    intersect = eye_o + d_v * t_c;
    if (!(intersect[0]*intersect[0] + intersect[2]*intersect[2] <= 0.25)) {
        t_c = -1;
    }
    return mPos(t_s, t_c);
}

float MyGLCanvas::mPos(float x, float y) {
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

double MyGLCanvas::intersectCylinder (glm::vec3 eyePointP, glm::vec3 rayV, glm::mat4 transformMatrix, glm::vec3 spherepos) {
	glm::mat4 transformInv = glm::inverse(transformMatrix);
	glm::vec3 eyePointObject = transformInv * glm::vec4(eyePointP, 1);
    glm::vec3 rayVObject = transformInv * glm::vec4(rayV, 0);
	glm::vec3 d = rayVObject;
	glm::vec3 eye = eyePointObject;
    double A = d[0]*d[0] + d[2]*d[2];
    double B = 2*eye[0]*d[0] + 2*eye[2]*d[2];
    double C = eye[0]*eye[0] + eye[2]*eye[2] - 0.25;
    float tb = quadraticForm(A, B, C);


    glm::vec3 intersect = eye + d * tb;
    if (!(intersect[1] > -0.5 && intersect[1] < 0.5)) {
        tb = -1;
    }

    float cap_one = (0.5 - eye[1]) / d[1];
	float cap_two = (-0.5 - eye[1]) / d[1];
    intersect = eye + d * cap_one;
    if (!(intersect[0]*intersect[0] + intersect[2]*intersect[2] <= 0.25)) {
        cap_one = -1;
    }
	intersect = eye + d * cap_two;
	if (!(intersect[0]*intersect[0] + intersect[2]*intersect[2] <= 0.25)) {
        cap_two = -1;
    }
    return mPos(tb, mPos(cap_one, cap_two));
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

SceneColor MyGLCanvas::computeColor(SceneMaterial material, glm::vec3 Nhat, glm::vec3 pos, glm::vec3 ray) {
	SceneGlobalData data;
	parser->getGlobalData(data);
	SceneColor color;
	color.r = 0;
	color.g = 0;
	color.b = 0;

	ray = glm::normalize(ray);

	// constant coeffecients we need
	float ka = data.ka;
	float ks = data.ks;
	float kd = data.kd;
	SceneColor Oa = material.cAmbient;
	SceneColor Os = material.cSpecular;
	SceneColor Od = material.cDiffuse;
	float shininess = material.shininess;
	// best attempt at doing the lighting equation
	// iterate through the lights and perform the equation in the handout
    for (int m = 0; m < parser->getNumLights(); m++) {
        SceneLightData lData;
        parser->getLightData(m, lData);
        SceneColor li = lData.color;
        glm::vec3 Lhati = lData.pos - pos;
        Lhati = glm::normalize(Lhati);
		glm::vec3 Ri = glm::normalize(Lhati - (2.0f * glm::dot(Lhati, Nhat) * Nhat));
		float RiV = glm::dot(Ri, glm::normalize(camera->getLookVector()));

        if (dot(Nhat, Lhati) > 0) {
			/*
            color.r += (kd * Od.r * li.r * dot(Nhat, Lhati));  
            color.g += (kd * Od.g * li.g * dot(Nhat, Lhati));  
            color.b += (kd * Od.b * li.b * dot(Nhat, Lhati));  
			color.r += (ks * Os.r) * pow(RiV, shininess);
			color.g += (ks * Os.g) * pow(RiV, shininess);
			color.b += (ks * Os.b) * pow(RiV, shininess);
			*/

			color.r += li.r * (kd * Od.r * dot(Nhat, Lhati) + (ks * Os.r) * pow(RiV, shininess));
			color.g += li.g * (kd * Od.g * dot(Nhat, Lhati) + (ks * Os.g) * pow(RiV, shininess));
			color.b += li.b * (kd * Od.b * dot(Nhat, Lhati) + (ks * Os.b) * pow(RiV, shininess));
        }
    }
    color.r += ka * (Oa.r);
    color.g += ka * (Oa.g);
    color.b += ka * (Oa.b);

	if (color.r * 255 > 255){
		color.r = 255;
	}
	else{
		color.r *= 255;
	}
    
    if (color.g * 255 > 255){
		color.g = 255;
	}
	else{
		color.g *= 255;
	}

	if (color.b * 255 > 255){
		color.b = 255;
	}
	else{
		color.b *= 255;
	}

	return color;
}

glm::vec3 MyGLCanvas::computeNormal(glm::vec3 inst, OBJ_TYPE shape) {
		if(shape == SHAPE_CYLINDER){
		    if (IN_RANGE(inst[1], 0.5))
               return glm::vec3(0, 1, 0);
           if (IN_RANGE(inst[1], -0.5))
               return glm::vec3(0, -1, 0);
			return glm::vec3(inst[0], 0, inst[2]);
		}else if (shape == SHAPE_CUBE){
			float inr = .4999;
            if (inst[0] > inr)
                return glm::vec3(1, 0, 0);
            if (inst[0] < -inr)
                return glm::vec3(-1, 0, 0);
            if (inst[1] > inr)
                return glm::vec3(0, 1, 0);
            if (inst[1] < -inr)
                return glm::vec3(0, -1, 0);
            if (inst[2] > inr)
                return glm::vec3(0, 0, 1);
            if (inst[2] < -inr)
                return glm::vec3(0, 0, -1);
		}else if (shape == SHAPE_SPHERE){
           return glm::vec3(inst[0], inst[1], inst[2]);
		}else if(shape == SHAPE_CONE){
			// base of cone
           if (IN_RANGE(inst[1], -0.5)){
               return glm::vec3(0, -1, 0);
		   }
           glm::vec3 vec1 = glm::vec3(inst[0], 0, inst[2]);
           vec1 = glm::normalize(vec1);
           glm::vec3 vec2 = glm::vec3(0, .5, 0);
           return vec1 + vec2; 
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
                bool hit = false;
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
								m = m * glm::translate(glm::mat4(1.0), my_trans->translate);
								break;
							case(TRANSFORMATION_SCALE):
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
						t = intersectCube(eye_pnt, ray, m, center);
						break;
					case SHAPE_CYLINDER:
						t = intersectCylinder(eye_pnt, ray, m, center);
						break;
					case SHAPE_CONE:
						t = intersectCone(eye_pnt, ray, m, center);
						break;
					case SHAPE_SPHERE:
						t = intersectSphere(eye_pnt, ray, m, center);
						break;
					}
				
					if (t > 0 && t < t_min) {
                        hit = true;

						t_min = t;
						intersection_obj = getIsectPointWorldCoord(glm::vec3(glm::inverse(m) * glm::vec4(eye_pnt,1.0)), glm::vec3(glm::inverse(m) * glm::vec4(ray,0)), t);
						glm::vec3 normal = computeNormal(intersection_obj, prim->type); 
						// normalize the normal
						normal = glm::normalize(glm::vec3(glm::transpose(glm::inverse(m)) * glm::vec4(normal,1)));
						glm::vec3 intersection = glm::vec3(transpose(m) * glm::vec4(intersection_obj, 1));
						color = computeColor(prim->material, normal, intersection, ray);
					}

				}

                if (hit) {
                    if (isectOnly == 1) {
                        setpixel(pixels, i, j, 255, 255, 255);
                    }
                    else {
                        setpixel(pixels, i, j, color.r, color.g, color.b);
                    }
                }
                else {
                    setpixel(pixels, i, j, 0, 0, 0);
                }

				// if (isectOnly == 1) {
				// 	setpixel(pixels, i, j, color.r, color.g, color.b);
				// }
				// else {
				// 	setpixel(pixels, i, j, 255, 255, 255);
				// }
			}
		}
	std::cout << "render complete" << endl;
	redraw();
}

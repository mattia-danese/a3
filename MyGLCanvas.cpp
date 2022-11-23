#define NUM_OPENGL_LIGHTS 8

#include "MyGLCanvas.h"
#include <glm/gtx/string_cast.hpp>
#include <float.h>
#include "random"

int Shape::m_segmentsX;
int Shape::m_segmentsY;

struct MyGLCanvas::nearestObj near_ist;

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

glm::vec3 MyGLCanvas::camLookAt(int pixelX, int pixelY){
	return glm::vec3(-1.0f + 2.0f * ((float)pixelX / (float)camera->getScreenWidth()),
		-1.0f + 2.0f * ((float)pixelY / (float)camera->getScreenHeight()),
		-1.0f);
}

glm::vec3 MyGLCanvas::generateRay(int pixelX, int pixelY) {
	lookAt = camLookAt(pixelX, pixelY);
	//likely stems from here
	glm::mat4 inv =  camera->getInverseModelViewMatrix() * camera->getInverseScaleMatrix();
	glm::vec4 worldFPP = inv * glm::vec4(lookAt, 1.0f);
	// already in world coordinates
	glm::vec4 worldEye = glm::vec4(getEyePoint(), 1.0f);

	glm::vec4 dHat = glm::normalize(worldFPP - worldEye);
	return glm::vec3(dHat.x, dHat.y, dHat.z);
    // lab has -dHat.y .. not sure which is best
}

glm::vec3 MyGLCanvas::getIsectPointWorldCoord(glm::vec3 eye, glm::vec3 ray, float t) {
	glm::vec3 p = eye + (t * ray);
	return p;
}

// FROM PREVIOUS LAB
double MyGLCanvas::intersectSphere (glm::vec3 eyePointP, glm::vec3 rayV, glm::mat4 transformMatrix, glm::vec3 spherepos) {
	glm::mat4 transformInv = glm::inverse(transformMatrix);
	glm::vec3 d = glm::vec3(transformInv * glm::vec4(rayV,0));
	glm::vec3 eye = glm::vec3(transformInv * glm::vec4(eyePointP, 1.0f));

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

float MyGLCanvas::quadraticForm(float A, float B, float C) {
    float det = B*B - 4.0*A*C;
	if (det < 0) {return -1;}
    else {
        double t1 = (-B + sqrt(det)) / (2.0 * A);
        double t2 = (-B - sqrt(det)) / (2.0 * A);

        if (t1 < 0) {return t2;}
        if (t2 < 0) {return t1;}
        return std::min(t1, t2);
    }
}

float MyGLCanvas::intersectCube (glm::vec3 eyePointP, glm::vec3 rayV, glm::mat4 transformMatrix, glm::vec3 spherepos) {
	glm::mat4 transformInv = glm::inverse(transformMatrix);
	glm::vec4 d_v = transformInv * glm::vec4(rayV,0);
	glm::vec4 eye = transformInv * glm::vec4(eyePointP, 1.0f);

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
	glm::vec3 eyePointObject = transformInv * glm::vec4(eyePointP, 1.0f);
    glm::vec3 rayVObject = transformInv * glm::vec4(rayV, 0);
	glm::vec3 d_v = rayVObject;
	glm::vec3 eye_o = eyePointObject;
    float A = d_v[0]*d_v[0] + d_v[2]*d_v[2] - (.25f*d_v[1]*d_v[1]);
    float B = 2.0f*eye_o[0]*d_v[0] + 2.0f*eye_o[2]*d_v[2] - .5f*eye_o[1]*d_v[1] + .25f*d_v[1];
    float C = eye_o[0]*eye_o[0] + eye_o[2]*eye_o[2] - .25f*pow(-eye_o[1]+.5, 2.0f);

    float t_s = quadraticForm(A, B, C);
    glm::vec3 intersect = eye_o + (d_v *t_s);
	if (intersect[1] < -0.5f || intersect[1] > 0.5f) {
        t_s = -1;
    }
	float t_c = (-0.5 - eye_o[1]) / d_v[1];
    glm::vec3 intersect_two = eye_o + (d_v * t_c);
    if (intersect_two[0]*intersect_two[0] + intersect_two[2]*intersect_two[2] > 0.25) {
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
	glm::vec3 eyePointObject = transformInv * glm::vec4(eyePointP, 1.0f);
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

SceneColor MyGLCanvas::bound(SceneColor c) {
    if (c.r > 1.0) c.r = 1.0;
    if (c.g > 1.0) c.g = 1.0;
    if (c.b > 1.0) c.b = 1.0;
	if (c.r < 0) c.r = 0.0;
	if (c.g < 0) c.g = 0.0;
	if (c.b < 0) c.b = 0.0;
    return c;
}

//vector<pair<ScenePrimitive*, vector<SceneTransformation*>>> my_scene_vals, glm::vec3 eye_pnt, glm::vec3 ray, int* hit
bool MyGLCanvas::shadowCheck(glm::vec3 pos, glm::vec3 Lhati){
	int hit = 0;
	// std::cout << "here" << std::endl;
	loopObjects(my_scene_vals, pos, Lhati, &hit, true);
	return hit;
}

SceneColor MyGLCanvas::computeColor(ScenePrimitive* p_m, glm::vec3 Nhat, glm::vec3 pos) {
	SceneMaterial material = p_m->material;
	SceneGlobalData data;
	parser->getGlobalData(data);
	SceneColor color;
	color.r = 0;
	color.g = 0;
	color.b = 0;

	// ray = glm::normalize(ray);

	glm::vec3 look = cameraData.look;
    look  = glm::normalize(look);
    glm::vec3 Vhat = pos - cameraData.pos;
    Vhat = glm::normalize(Vhat);

	// constant coeffecients we need
	float ka = data.ka;
	float ks = data.ks;
	float kd = data.kd;
	SceneColor Oa = bound(material.cAmbient);
	SceneColor Os = bound(material.cSpecular);
	SceneColor Od = bound(material.cDiffuse);
	SceneColor Or = bound(material.cReflective);
	float shininess = material.shininess;
	// best attempt at doing the lighting equation
	// iterate through the lights and perform the equation in the handout
	SceneLightData lData;

		SceneColor Ir;
    if (depth_fresh > 0 && !material.textureMap->isUsed){
		// std::cout << "recursing" << std::endl;
        depth_fresh--;
        glm::vec3 d = Vhat - 2 * (glm::dot(Vhat, Nhat)) * Nhat;
        d = glm::normalize(d);
		int hit = 0;
		//Ir = findColor(findNearestIntersection(pos, d), pos, d);
        Ir = loopObjects(my_scene_vals, pos, d, &hit, false);
		Ir = bound(Ir);
    } else { 
        Ir.r = 0;
        Ir.g = 0;
        Ir.b = 0;
        Ir.a = 0;
    }

    for (int m = 0; m < parser->getNumLights(); m++) {
		parser->getLightData(m, lData);
        SceneColor li = bound(lData.color);
        glm::vec3 Lhati;
		if(lData.type == LIGHT_DIRECTIONAL){
			Lhati = glm::vec3(-lData.dir.x, -lData.dir.y, -lData.dir.z);
		}
		if(lData.type == LIGHT_POINT){
			Lhati = lData.pos - pos;
		}
		if(shadowCheck(pos, Lhati)){
			continue;
		}
        Lhati = glm::normalize(Lhati);
		glm::vec3 Lvec  = pos - lData.pos;
		glm::vec3 Ri = glm::normalize(Lhati - (2.0f * glm::dot(Lhati, Nhat) * Nhat));
		float RiV = glm::dot(Ri, glm::normalize(camera->getLookVector()));
        if (dot(Nhat, Lhati) > 0) {
			color.r += li.r * (kd * Od.r * dot(Nhat, Lhati) + (ks * Os.r) * pow(RiV, shininess));
			color.g += li.g * (kd * Od.g * dot(Nhat, Lhati) + (ks * Os.g) * pow(RiV, shininess));
			color.b += li.b * (kd * Od.b * dot(Nhat, Lhati) + (ks * Os.b) * pow(RiV, shininess));
            }
    }
	color = bound(color);

	// color = textureMap(color, prim_m);

    color.r += ka * (Oa.r) + ks * Or.r * Ir.r;
    color.g += ka * (Oa.g)+ ks * Or.g * Ir.g;
    color.b += ka * (Oa.b)+  ks * Or.b * Ir.b;
	color = bound(color);

	// color.r += ka * (Oa.r);
    // color.g += ka * (Oa.g);
    // color.b += ka * (Oa.b);

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


void MyGLCanvas::convert_xyz_to_cube_uv(float x, float y, float z, int *index, float *u, float *v)
{
  float absX = fabs(x);
  float absY = fabs(y);
  float absZ = fabs(z);
  
  int isXPositive = x > 0 ? 1 : 0;
  int isYPositive = y > 0 ? 1 : 0;
  int isZPositive = z > 0 ? 1 : 0;
  
  float maxAxis, uc, vc;
  
  // POSITIVE X
  if (isXPositive && absX >= absY && absX >= absZ) {
    // u (0 to 1) goes from +z to -z
    // v (0 to 1) goes from -y to +y
    maxAxis = absX;
    uc = -z;
    vc = y;
    *index = 0;
  }
  // NEGATIVE X
  if (!isXPositive && absX >= absY && absX >= absZ) {
    // u (0 to 1) goes from -z to +z
    // v (0 to 1) goes from -y to +y
    maxAxis = absX;
    uc = z;
    vc = y;
    *index = 1;
  }
  // POSITIVE Y
  if (isYPositive && absY >= absX && absY >= absZ) {
    // u (0 to 1) goes from -x to +x
    // v (0 to 1) goes from +z to -z
    maxAxis = absY;
    uc = x;
    vc = z;
    *index = 2;
  }
  // NEGATIVE Y
  if (!isYPositive && absY >= absX && absY >= absZ) {
    // u (0 to 1) goes from -x to +x
    // v (0 to 1) goes from -z to +z
    maxAxis = absY;
    uc = x;
    vc = -z;
    *index = 3;
  }
  // POSITIVE Z
  if (isZPositive && absZ >= absX && absZ >= absY) {
    // u (0 to 1) goes from -x to +x
    // v (0 to 1) goes from -y to +y
    maxAxis = absZ;
    uc = x;
    vc = y;
    *index = 4;
  }
  // NEGATIVE Z
  if (!isZPositive && absZ >= absX && absZ >= absY) {
    // u (0 to 1) goes from +x to -x
    // v (0 to 1) goes from -y to +y
    maxAxis = absZ;
    uc = -x;
    vc = y;
    *index = 5;
  }

  // Convert range from -1 to 1 to 0 to 1
  *u = .5f * (uc / maxAxis + 1.0f);
  *v = .5f * (vc / maxAxis + 1.0f);
}



SceneColor MyGLCanvas::textureMap(SceneColor color, ScenePrimitive* prim){
		float blend = 0;
        glm::vec3 color_n;
        SceneColor tex_color_blend;
        if(prim != nullptr && prim->material.textureMap->isUsed && !p_map.empty() && ist_min != glm::vec3(INFINITY)){
                    SceneMaterial mat = prim->material;
                    SceneFileMap* texture = mat.textureMap;
                    ppm* my_ppm;
                    my_ppm = p_map[prim_m->material.textureMap->filename];
					// for(int i = 0; i < my_ppm->getWidth(); i++){
					// 	for(int j = 0; j < my_ppm->getHeight(); j++){
					// 		SceneColor px = my_ppm->getPixel(i, my_ppm->getHeight()-j);
					// 		setpixel(pixels, i, j, px.r, px.g, px.b);
					// 		//setpixel(pixels, i,  my_ppm->getHeight()-j, px.r, px.g, px.b);
					// 	}
					// }
                    blend = mat.blend;
                    float u = 0;
                    float v = 0;
                    float s = 0;
                    float t = 0;
                    glm::vec3 v_cor;
                    glm::vec3 u_cor;
                    glm::vec3 n;
                    float rawU;
                    float theta;
					float phi;
                    int index = 0;
                    switch (prim->type) {
                    case SHAPE_CUBE:
                        convert_xyz_to_cube_uv(ist_min.x, ist_min.y, ist_min.z, &index, &u, &v);
                        break;
                    case SHAPE_CYLINDER:
                        theta = .5 + atan2(ist_min.x, ist_min.z);
                        rawU = theta / (2.0f * PI);
                        u = my_ppm->getWidth()-(1 - (rawU + 0.5));
                        v = fmod(.5 + ist_min[1], 1);
                        //cone cap
                        if(ist_min.y > .499 || ist_min.y < -.499){
							// u = atan2(ist_min.y, ist_min.x);
                            // if(theta < 0){
                            //     u = -theta / (2.0f * PI);
                            // }else{
                            //     u = theta / (2.0f * PI);
                            // }
							// float radius = sqrt(pow(ist_min.x, 2)+ pow(ist_min.y, 2)+ pow(ist_min.z, 2));
							// phi = acos(ist_min.y / radius);
							// v = 1 - (phi/PI);
                            convert_xyz_to_cube_uv(ist_min.x, ist_min.y, ist_min.z, &index, &u, &v);
                        }
                        break;
                    case SHAPE_CONE:
					    theta = .5 + atan2(ist_min.x, ist_min.z);
                        rawU = theta / (2.0f * PI);
                        u = 1 - (rawU + 0.5);
                        v = fmod(.5 + ist_min[1], 1);
                        //cone cap
                        if(ist_min.y > .499 || ist_min.y < -.499){
                            // if(theta < 0){
                            //     u = -theta / (2.0f * PI);
                            // }else{
                            //     u = theta / (2.0f * PI);
                            // }
							// float radius = sqrt(pow(ist_min.x, 2)+ pow(ist_min.y, 2)+ pow(ist_min.z, 2));
							// phi = acos(ist_min.y / radius);
							// v = 1 - (phi/PI);
                            convert_xyz_to_cube_uv(ist_min.x, ist_min.y, ist_min.z, &index, &u, &v);
                        }
                        break;
                    case SHAPE_SPHERE:
						float radius = sqrt(pow(ist_min.x, 2)+ pow(ist_min.y, 2)+ pow(ist_min.z, 2));
						phi = acos(ist_min.y / radius);
                        u = 1 - (atan2(ist_min.x, ist_min.z) / (2.0f * PI));
                        v = 1 - (phi/PI);
                        break;
                    }
					
                    s = fmod((u*(float)my_ppm->getWidth()*texture->repeatU), (float)my_ppm->getWidth());
                    t = fmod((v*(float)my_ppm->getHeight()*texture->repeatV), (float)my_ppm->getHeight());
					SceneColor tex_c = my_ppm->getPixel(s, my_ppm->getHeight()-t);
                    if(tex_c.r == 0 && tex_c.g == 0 && tex_c.b == 0){
                        return color;
                    }
                    color_n = glm::vec3(color.r, color.g, color.b) * (1-blend) + glm::vec3(tex_c.r, tex_c.g, tex_c.b)*blend;
                }else{
                    return color;
                }
                tex_color_blend.r = color_n.x;
                tex_color_blend.g = color_n.y;
                tex_color_blend.b = color_n.z;
                return tex_color_blend;
}

glm::vec3 MyGLCanvas::computeNormal(glm::vec3 inst, OBJ_TYPE shape) {
	// just a small floating point value to capture the base of the objects
	// in object coordinate space
	float inr = .4999;
		if(shape == SHAPE_CYLINDER){
            if (inst[1] > inr)
                return glm::vec3(0, 1, 0);
            if (inst[1] < -inr)
                return glm::vec3(0, -1, 0);
			return glm::vec3(inst[0], 0, inst[2]);
		}else if (shape == SHAPE_CUBE){
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
		    if (inst[1] < -inr)
                 return glm::vec3(0, -1, 0);
           glm::vec3 vec1 = glm::vec3(inst[0], 0, inst[2]);
           vec1 = glm::normalize(vec1);
           glm::vec3 vec2 = glm::vec3(0, .5, 0);
           return vec1 + vec2; 
		   }
}



void MyGLCanvas::traverse1(SceneNode* root, vector<pair<ScenePrimitive*, vector<SceneTransformation*>>>& my_scene_vals, vector<SceneTransformation*> curr_trans, map<string, ppm*>* p_m)
{
	if (root == nullptr)
		return;

	for (SceneNode* node : root->children) { // storing transformations in non-primitive nodes 
		for (SceneTransformation* my_trans : node->transformations) {
			curr_trans.push_back(my_trans);
		}

		traverse1(node, my_scene_vals, curr_trans, p_m);
		for (ScenePrimitive* my_prim : node->primitives) { // for all primitives in leaf node, make pair of prim and previously collected trasnformations
			list<SceneTransformation*> curr_scene_trans = {};
			if(my_prim->material.textureMap->isUsed){
				std::cout << "here?" << std::endl;
				std::cout << my_prim->material.textureMap->isUsed << std::endl;
				std::cout << my_prim->material.textureMap->filename << std::endl;
				if (p_m->find(my_prim->material.textureMap->filename) == p_m->end()) {
					p_m->insert(pair<string, ppm*>(my_prim->material.textureMap->filename, new ppm(my_prim->material.textureMap->filename)));
				}
			}
			my_scene_vals.push_back(pair<ScenePrimitive*, vector<SceneTransformation*>>(my_prim, curr_trans));
		}
		for (SceneTransformation* my_trans : node->transformations) { // func will recurse back up, delete transformations
			curr_trans.pop_back();
		}
	}
}

SceneColor MyGLCanvas::loopObjects(vector<pair<ScenePrimitive*, vector<SceneTransformation*>>> my_scene_vals, glm::vec3 eye_pnt, glm::vec3 ray, int* hit, bool shadow_check){
	        glm::vec3 intersection_obj = glm::vec3(0);
			glm::vec3 intersection = glm::vec3(0);
			glm::vec3 normal = glm::vec3(0);
			SceneColor color;

			ScenePrimitive* prim = nullptr;
			color.r = 0, color.g = 0, color.b = 0;
			*hit = false;
			float t_min = FLT_MAX;
			for (pair<ScenePrimitive*, vector<SceneTransformation*>> my_p : my_scene_vals) {
				// std::cout << "here" << std::endl;
					prim = my_p.first;
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
								//std::cout << my_trans->angle << std::endl;
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
				
					if (t > 0.01 && (t_min < 0 || t < t_min)) {
                        *hit = true;
						if(shadow_check){
							return color;
						}
						t_min = t;
						prim_m = prim;
						//object coordinates intersection
						intersection_obj = getIsectPointWorldCoord(glm::vec3(glm::inverse(m) * glm::vec4(eye_pnt,1.0f)), glm::vec3(glm::inverse(m) * glm::vec4(ray,0)), t_min);
						ist_min = intersection_obj; 
						normal = computeNormal(intersection_obj, prim->type);
						normal = glm::normalize(glm::vec3(glm::transpose(glm::inverse(m)) * glm::vec4(normal,1))); 
						intersection = glm::vec3(m * glm::vec4(intersection_obj, 1));
					}
				}
				if(t_min != FLT_MAX && prim != nullptr){
					//std::cout << "here" << std::endl;
					color = textureMap(computeColor(prim_m, normal, intersection), prim_m);
				}
				return color;
}


inline double random_double() {
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return distribution(generator);
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
	ist_min = glm::vec3(INFINITY);

	
	pixels = new GLubyte[pixelWidth  * pixelHeight * 3];
	memset(pixels, 0, pixelWidth  * pixelHeight * 3);
		glm::vec3  eye_pnt = getEyePoint();
				if (my_scene_vals.empty()) {
					std::cout << "calling traverse!" << endl;
					SceneNode* root = parser->getRootNode();
					traverse1(root, my_scene_vals, scenetransformations, &p_map);
				}
			std::cout << "traverse done!" << endl;
		for (int i = 0; i < pixelWidth; i++) {
			for (int j = 0; j < pixelHeight; j++) {
				int hit = 0;
				SceneColor color;
				color.r = 0;
				color.g = 0;
				color.b = 0; 
				SceneColor temp;
				for (int s = 0; s < samples_per_pixel; ++s) {
				float t_min = FLT_MAX;
				color.r = 0, color.g = 0, color.b = 0;
				//TODO: this is where your ray casting will happen!
				glm::vec3 ray = generateRay(i, j);
				depth_fresh = depth;
				// color = findColor(findNearestIntersection(eye_pnt, ray), eye_pnt, ray);
				// setpixel(pixels, i, j, color.r, color.g, color.b);

				temp = loopObjects(my_scene_vals, eye_pnt, ray, &hit, false);
				color.r += temp.r;
				color.g += temp.g;
				color.b += temp.b;
				}
				float scale = 1.0f/(float)samples_per_pixel;
				color.r *= scale;
				color.g *= scale;
				color.b *= scale;
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
			}
		}
	std::cout << "render complete" << endl;
	redraw();
}

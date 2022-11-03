#include "Camera.h"
#include "glm/glm.hpp"
#include <glm/gtc/type_ptr.hpp>

Camera::Camera() {
	reset();
}

Camera::~Camera() {
}

glm::vec3 Camera::cross_product_n(glm::vec3 first, glm::vec3 second) {
	// find vector between x2 and x1
	float v1x = first[0];
	float v1y = first[1];
	float v1z = first[2];

	// find vector between x3 and x2
	float v2x = second[0];
	float v2y = second[1];
	float v2z = second[2];

	// cross product v1xv2

	float cx = v1y * v2z - v1z * v2y;
	float cy = v1z * v2x - v1x * v2z;
	float cz = v1x * v2y - v1y * v2x;

	// normalize

	float length = sqrt(cx * cx + cy * cy + cz * cz);
	if (length != 0) {
		cx = cx / length;
		cy = cy / length;
		cz = cz / length;
	}
	return glm::vec3(cx, cy, cz);
}

glm::vec3 Camera::cross_product(glm::vec3 first, glm::vec3 second) {
	// find vector between x2 and x1
	float v1x = first[0];
	float v1y = first[1];
	float v1z = first[2];

	// find vector between x3 and x2
	float v2x = second[0];
	float v2y = second[1];
	float v2z = second[2];

	// cross product v1xv2

	float cx = v1y * v2z - v1z * v2y;
	float cy = v1z * v2x - v1x * v2z;
	float cz = v1x * v2y - v1y * v2x;

	return glm::vec3(cx, cy, cz);
}


void Camera::reset() {
	setViewAngle(VIEW_ANGLE);
	setNearPlane(NEAR_PLANE);
	setFarPlane(FAR_PLANE);
	screenWidth = screenHeight = 200;
	screenWidthRatio = 1.0f;
	rotU = rotV = rotW = 0;
	orientLookAt(glm::vec3(0.0f, 0.0f, DEFAULT_FOCUS_LENGTH), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));


}

//called by main.cpp as a part of the slider callback for controlling rotation
// the reason for computing the diff is to make sure that we are only incrementally rotating the camera
void Camera::setRotUVW(float u, float v, float w) {
	float diffU = u - rotU;
	float diffV = v - rotV;
	float diffW = w - rotW;
	rotateU(diffU);
	rotateV(diffV);
	rotateW(diffW);
	rotU = u;
	rotV = v;
	rotW = w;
}


void Camera::orientLookAt(glm::vec3 eyePoint, glm::vec3 lookatPoint, glm::vec3 upVec) {
	this->w = -normalize(lookatPoint - eyePoint);
	this->u = normalize(cross(upVec, w));
	this->v = cross(this->w, this->u);

	this->eye = eyePoint;
	//wuv maps to fsu
	this->lookVec = lookVec;
	this->upVec = upVec;
	calculate_matrices();
}


void Camera::orientLookVec(glm::vec3 eyePoint, glm::vec3 lookVec, glm::vec3 upVec) {
	this->w = -normalize(lookVec - eyePoint);
	this->u = normalize(cross(upVec, w));
	this->v = cross(this->w, this->u);
	this->eye = eyePoint;
	//wuv maps to fsu
	this->lookVec = lookVec;
	this->upVec = upVec;
	calculate_matrices();

}

void Camera::calculate_matrices() {
	this->rotate_inverse = {
		this->u[0], this->v[0], this->w[0], 0.0f,
		this->u[1], this->v[1], this->w[1], 0.0f,
		this->u[2], this->v[2], this->w[2], 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	};

	this->translate_inverse = {
		1.0f, 0, 0, 0,
		0, 1.0f, 0, 0,
		0, 0, 1.0f, 0,
		-this->eye[0], -this->eye[1], -this->eye[2], 1.0f
	};

	float theta_w = glm::radians(this->viewAngle);
	float h2 = tan(theta_w / 2.0f) * farPlane * (1.0f + this->getScreenHeight()) / (1.0f + this->getScreenWidth());
	float theta_h = atan(h2 / this->farPlane) * 2.0f;
	float w = 1.0f / (tan(theta_w / 2.0f) * this->farPlane);
	float h = 1.0f / (tan(theta_h / 2.0f) * this->farPlane);
	this->scale = {
		w, 0, 0, 0,
		0, h, 0, 0,
		0, 0, 1.0f / this->farPlane, 0,
		0, 0, 0, 1.0f
	};


	float c = -this->nearPlane / this->farPlane;
	this->m_pp = {
		1.0f, 0, 0, 0,
		0, 1.0f, 0, 0,
		0, 0, -1.0f / (c + 1.0f), -1.0f,
		0, 0,  c / (c + 1.0f), 0
	};
	this->model_view_matrix = multiply_matrix(this->translate_inverse, this->rotate_inverse);
	this->projection_matrix = multiply_matrix(this->scale, this->m_pp);
	float end = 0;
}

glm::vec4 Camera::multiply_vector(glm::mat4 m, glm::vec4 p) {
	glm::vec4 v_ret;
	v_ret[0] = p[0] * m[0][0] + p[1] * m[0][1] + p[2] * m[0][2] + p[3] * m[0][3];
	v_ret[1] = p[0] * m[1][0] + p[1] * m[1][1] + p[2] * m[1][2] + p[3] * m[1][3];
	v_ret[2] = p[0] * m[2][0] + p[1] * m[2][1] + p[2] * m[2][2] + p[3] * m[2][3];
	v_ret[3] = p[0] * m[3][0] + p[1] * m[3][1] + p[2] * m[3][2] + p[3] * m[3][3];
	return v_ret;

}

glm::mat4 Camera::multiply_matrix(glm::mat4 m1, glm::mat4 m2) {
	return glm::mat4(
		m1[0][0] * m2[0][0] + m1[0][1] * m2[1][0] + m1[0][2] * m2[2][0] + m1[0][ 3] * m2[3][ 0],
		m1[0][ 0] * m2[0][ 1] + m1[0][ 1] * m2[1][ 1] + m1[0][ 2] * m2[2][ 1] + m1[0][ 3] * m2[3][ 1],
		m1[0][ 0] * m2[0][ 2] + m1[0][ 1] * m2[1][ 2] + m1[0][ 2] * m2[2][ 2] + m1[0][ 3] * m2[3][ 2],
		m1[0][ 0] * m2[0][ 3] + m1[0][ 1] * m2[1][ 3] + m1[0][ 2] * m2[2][ 3] + m1[0][ 3] * m2[3][ 3],

		m1[1][ 0] * m2[0][ 0] + m1[1][ 1] * m2[1][ 0] + m1[1][ 2] * m2[2][ 0] + m1[1][ 3] * m2[3][ 0],
		m1[1][ 0] * m2[0][ 1] + m1[1][ 1] * m2[1][ 1] + m1[1][ 2] * m2[2][ 1] + m1[1][ 3] * m2[3][ 1],
		m1[1][ 0] * m2[0][ 2] + m1[1][ 1] * m2[1][ 2] + m1[1][ 2] * m2[2][ 2] + m1[1][ 3] * m2[3][ 2],
		m1[1][ 0] * m2[0][ 3] + m1[1][ 1] * m2[1][ 3] + m1[1][ 2] * m2[2][ 3] + m1[1][ 3] * m2[3][ 3],

		m1[2][ 0] * m2[0][ 0] + m1[2][ 1] * m2[1][ 0] + m1[2][ 2] * m2[2][ 0] + m1[2][ 3] * m2[3][ 0],
		m1[2][ 0] * m2[0][ 1] + m1[2][ 1] * m2[1][ 1] + m1[2][ 2] * m2[2][ 1] + m1[2][ 3] * m2[3][ 1],
		m1[2][ 0] * m2[0][ 2] + m1[2][ 1] * m2[1][ 2] + m1[2][ 2] * m2[2][ 2] + m1[2][ 3] * m2[3][ 2],
		m1[2][ 0] * m2[0][ 3] + m1[2][ 1] * m2[1][ 3] + m1[2][ 2] * m2[2][ 3] + m1[2][ 3] * m2[3][ 3],

		m1[3][ 0] * m2[0][ 0] + m1[3][ 1] * m2[1][ 0] + m1[3][ 2] * m2[2][ 0] + m1[3][ 3] * m2[3][ 0],
		m1[3][ 0] * m2[0][ 1] + m1[3][ 1] * m2[1][ 1] + m1[3][ 2] * m2[2][ 1] + m1[3][ 3] * m2[3][ 1],
		m1[3][ 0] * m2[0][ 2] + m1[3][ 1] * m2[1][ 2] + m1[3][ 2] * m2[2][ 2] + m1[3][ 3] * m2[3][ 2],
		m1[3][ 0] * m2[0][ 3] + m1[3][ 1] * m2[1][ 3] + m1[3][ 2] * m2[2][ 3] + m1[3][ 3] * m2[3][ 3]
	);
}



glm::mat4 Camera::getScaleMatrix() {
	return this->scale;
}

glm::mat4 Camera::getInverseScaleMatrix() {
	return glm::inverse(this->scale);
}

glm::mat4 Camera::getUnhingeMatrix() {
	return this->m_pp;
}


glm::mat4 Camera::getProjectionMatrix() {
	return this->projection_matrix;
}

glm::mat4 Camera::getInverseModelViewMatrix() {
	return this->rotate * this->translate;
}


void Camera::setViewAngle (float _viewAngle) {
	this->viewAngle = _viewAngle;
}

void Camera::setNearPlane (float _nearPlane) {
	this->nearPlane = _nearPlane;
}

void Camera::setFarPlane (float _farPlane) {
	this->farPlane = _farPlane;
}

void Camera::setScreenSize (int _screenWidth, int _screenHeight) {
	this->screenWidth = _screenWidth;
	this->screenHeight = _screenHeight;
	this->screenWidthRatio = _screenWidth / _screenHeight;
}

glm::mat4 Camera::getModelViewMatrix() {
	return this->model_view_matrix;
}


void Camera::rotateV(float degrees) {
	glm::mat4 mymatrix = glm::rotate(glm::mat4(1.0), -glm::radians(degrees), this->v);
	this->u = mymatrix * glm::vec4(this->u, 1);
	this->v = mymatrix * glm::vec4(this->v, 1);
	this->w = mymatrix * glm::vec4(this->w, 1);
	calculate_matrices();

}

void Camera::rotateU(float degrees) {
	glm::mat4 mymatrix = glm::rotate(glm::mat4(1.0), glm::radians(degrees), this->u);
	this->u = mymatrix * glm::vec4(this->u, 1);
	this->v = mymatrix * glm::vec4(this->v, 1);
	this->w = mymatrix * glm::vec4(this->w, 1);
	calculate_matrices();


}

void Camera::rotateW(float degrees) {
	glm::mat4 mymatrix = glm::rotate(glm::mat4(1.0), -glm::radians(degrees), this->w);
	this->u = mymatrix * glm::vec4(this->u, 1);
	this->v = mymatrix * glm::vec4(this->v, 1);
	this->w = mymatrix * glm::vec4(this->w, 1);
	calculate_matrices();
}

//void Camera::translate(glm::vec3 v) {
//}
//
//void Camera::rotate(glm::vec3 point, glm::vec3 axis, float degrees) {
//}


glm::vec3 Camera::getEyePoint() {
	return this->eye;
}

glm::vec3 Camera::getLookVector() {
	return -this->w;
}

glm::vec3 Camera::getUpVector() {
	return this->upVec;
}

float Camera::getViewAngle() {
	return viewAngle;
}

float Camera::getNearPlane() {
	return nearPlane;
}

float Camera::getFarPlane() {
	return farPlane;
}

int Camera::getScreenWidth() {
	return screenWidth;
}

int Camera::getScreenHeight() {
	return screenHeight;
}

float Camera::getScreenWidthRatio() {
	return this->screenWidth/this->screenHeight;
}

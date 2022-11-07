#include "Camera.h"
#include <glm/gtx/string_cast.hpp>
#include <iostream>

Camera::Camera() {
	reset();
}

Camera::~Camera() {
}

void Camera::reset() {
	orientLookAt(glm::vec3(0.0f, 0.0f, DEFAULT_FOCUS_LENGTH), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
	setViewAngle(VIEW_ANGLE);
	setNearPlane(NEAR_PLANE);
	setFarPlane(FAR_PLANE);
	screenWidth = screenHeight = 200;
	screenWidthRatio = 1.0f;
	rotU = rotV = rotW = 0;

	updateProject = true;
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
	eyePos = eyePoint;
	lookVector = glm::normalize(glm::vec3(lookatPoint.x - eyePoint.x, lookatPoint.y - eyePoint.y, lookatPoint.z - eyePoint.z));
	glm::vec3 rightVec = glm::normalize(glm::cross(lookVector, glm::vec3(0.0f, 1.0f, 0.0f)));
	upVector = glm::normalize(glm::cross(rightVec, lookVector));

	tranlateM = glm::translate(glm::mat4(1.0), glm::vec3(-1 * eyePoint.x, -1 * eyePoint.y, -1 * eyePoint.z));
	updateRotateMatrix(lookVector, upVector);
}


void Camera::orientLookVec(glm::vec3 eyePoint, glm::vec3 lookVec, glm::vec3 upVec) {
	eyePos = eyePoint;
	lookVector = glm::normalize(lookVec);
	// glm::vec3 rightVec = glm::normalize(glm::cross(lookVector, glm::vec3(0.0f, 1.0f, 0.0f)));
    // upVector = glm::normalize(glm::cross(rightVec, lookVector));
    upVector = upVec;
    // std::cout << glm::to_string(upVector) << std::endl;

	tranlateM = glm::translate(glm::mat4(1.0), glm::vec3(-1 * eyePoint.x, -1 * eyePoint.y, -1 * eyePoint.z));
	updateRotateMatrix(lookVector, upVector);
}

void Camera::updateRotateMatrix(glm::vec3 lookVec, glm::vec3 upVec) {
	w = glm::normalize(glm::vec3(-1 * lookVec.x, -1 * lookVec.y, -1 * lookVec.z));
	u = glm::normalize(glm::cross(upVec, w));
	v = glm::normalize(glm::cross(w, u));
	rotateM = glm::mat4(1.0);

	rotateM[0][0] = u.x;
	rotateM[1][0] = u.y;
	rotateM[2][0] = u.z;

	rotateM[0][1] = v.x;
	rotateM[1][1] = v.y;
	rotateM[2][1] = v.z;

	rotateM[0][2] = w.x;
	rotateM[1][2] = w.y;
	rotateM[2][2] = w.z;

	uvw = rotateM;
}

glm::mat4 Camera::getScaleMatrix() {
	if (updateProject)
		this->updateProjectMatrix();
	return scaleM;
}

glm::mat4 Camera::getInverseScaleMatrix() {
	if (updateProject)
		this->updateProjectMatrix();
	glm::mat4 temp = glm::mat4(1.0);
	temp[0][0] = 1 / scaleM[0][0];
	temp[1][1] = 1 / scaleM[1][1];
	temp[2][2] = 1 / scaleM[2][2];
	return temp;
}

glm::mat4 Camera::getUnhingeMatrix() {
	if (updateProject)
		this->updateProjectMatrix();
	return mppM;
}


glm::mat4 Camera::getProjectionMatrix() {
	if (updateProject)
		this->updateProjectMatrix();
	glm::mat4 projectM = mppM * scaleM;
	return projectM;
}

glm::mat4 Camera::getInverseModelViewMatrix() {
	glm::mat4 tInverse = glm::mat4(1.0);
	tInverse[3][0] = -1 * tranlateM[3][0];
	tInverse[3][1] = -1 * tranlateM[3][1];
	tInverse[3][2] = -1 * tranlateM[3][2];

	glm::mat4 rInverse = glm::transpose(rotateM);
	return tInverse * rInverse;
}

glm::mat4 Camera::getModelViewMatrix() {
	glm::mat4 modelViewMat4 = rotateM * tranlateM;
	return modelViewMat4;
}


void Camera::setViewAngle (float _viewAngle) {
	viewAngle = _viewAngle;
	updateProject = true;
}

void Camera::setNearPlane (float _nearPlane) {
	nearPlane = _nearPlane;
	filmPlanDepth = _nearPlane;
	updateProject = true;
}

void Camera::setFarPlane (float _farPlane) {
	farPlane = _farPlane;
	updateProject = true;
}

void Camera::setScreenSize (int _screenWidth, int _screenHeight) {
	screenWidth = _screenWidth;
	screenHeight = _screenHeight;
	screenWidthRatio = (float)_screenWidth / (float)_screenHeight;
	updateProject = true;
}

void Camera::updateProjectMatrix() {
	//scale matrix
	float widthAngle = viewAngle * screenWidthRatio;
	scaleM = glm::mat4(1.0f);
	scaleM[0][0] = 1.0 / (tan((widthAngle / 2) * PI / 180.0) * farPlane);
	scaleM[1][1] = 1.0 / (tan((viewAngle / 2) * PI / 180.0) * farPlane);
	scaleM[2][2] = 1.0 / farPlane;

	//unhinge matrix
	float c = -1 * (nearPlane / farPlane);
	mppM = glm::mat4(1.0f);
	mppM[2][2] = -1 / (c + 1);
	mppM[2][3] = -1;
	mppM[3][2] = c / (c + 1);
	mppM[3][3] = 0;

	updateProject = false;
}


void Camera::rotateV(float degrees) {
	rotate(glm::vec3(), v, degrees);
}

void Camera::rotateU(float degrees) {
	rotate(glm::vec3(), u, degrees);
}

void Camera::rotateW(float degrees) {
	rotate(glm::vec3(), glm::vec3(-1 * w.x, -1 * w.y, -1 * w.z), degrees);
}

void Camera::translate(glm::vec3 v) {
	tranlateM = glm::translate(glm::mat4(1.0), glm::vec3(-1 * v.x, -1 * v.y, -1 * v.z));
}

void Camera::rotate(glm::vec3 point, glm::vec3 axis, float degrees) {
	glm::mat4 trans = glm::mat4(1.0f);
	trans = glm::rotate(trans, glm::radians(degrees), axis);
	lookVector = trans * glm::vec4(lookVector, 1.0);
	upVector = trans * glm::vec4(upVector, 1.0);
	updateRotateMatrix(lookVector, upVector);
}


glm::vec3 Camera::getEyePoint() {
	return eyePos;
}

glm::vec3 Camera::getLookVector() {
	return lookVector;
}

glm::vec3 Camera::getUpVector() {
	return upVector;
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
	return screenWidthRatio;
}
#include "Camera.h"

Camera::Camera() {
	reset();
}

Camera::~Camera() {
}

void Camera::reset() {
    eyePoint = glm::vec3();
    basis = glm::mat4(1.0);
	orientLookAt(glm::vec3(0.0f, 0.0f, DEFAULT_FOCUS_LENGTH), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
	setViewAngle(VIEW_ANGLE);
	setNearPlane(NEAR_PLANE);
	setFarPlane(FAR_PLANE);
	screenWidth = screenHeight = 200;
	screenWidthRatio = 1.0f;
	rotU = rotV = rotW = 0;
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
    orientLookVec(eyePoint, lookatPoint - eyePoint, upVec);
}


void Camera::orientLookVec(glm::vec3 eyePoint, glm::vec3 lookVec, glm::vec3 upVec) {
    this->eyePoint = eyePoint;

    glm::vec3 w = -glm::normalize(lookVec);
    glm::vec3 u = glm::cross(upVec, w);
    u = glm::normalize(u);
    glm::vec3 v = glm::cross(w, u);
    
    basis[0][0] = u.x;
    basis[0][1] = u.y;
    basis[0][2] = u.z;
    basis[1][0] = v.x;
    basis[1][1] = v.y;
    basis[1][2] = v.z;
    basis[2][0] = w.x;
    basis[2][1] = w.y;
    basis[2][2] = w.z;
}

glm::mat4 Camera::getScaleMatrix() {
    float widthAngle = viewAngle * screenWidthRatio;
    float xScale = 1 / (glm::tan(glm::radians(widthAngle) / 2) * farPlane);
    float yScale = 1 / (glm::tan(glm::radians(viewAngle) / 2) * farPlane);
    float zScale = 1 / farPlane;
    return glm::scale(glm::mat4(1.0), glm::vec3(xScale, yScale, zScale));
}

glm::mat4 Camera::getInverseScaleMatrix() {
    float widthAngle = viewAngle * screenWidthRatio;
    float xScale = glm::tan(glm::radians(widthAngle) / 2) * farPlane;
    float yScale = glm::tan(glm::radians(viewAngle) / 2) * farPlane;
    float zScale = farPlane;
    return glm::scale(glm::mat4(1.0), glm::vec3(xScale, yScale, zScale));
}

glm::mat4 Camera::getUnhingeMatrix() {
    float c = -nearPlane / farPlane;
	glm::mat4 unhingeMat(1.0);
    unhingeMat[2][2] = -1 / (c + 1);
    unhingeMat[3][2] = c / (c + 1);
    unhingeMat[3][3] = 0;
    unhingeMat[2][3] = -1;
	return unhingeMat;
}


glm::mat4 Camera::getProjectionMatrix() {
	return getUnhingeMatrix() * getScaleMatrix();
}

glm::mat4 Camera::getInverseModelViewMatrix() {
    return glm::translate(glm::mat4(1.0), eyePoint) * basis;
}


void Camera::setViewAngle (float _viewAngle) {
    viewAngle = _viewAngle;
}

void Camera::setNearPlane (float _nearPlane) {
    nearPlane = _nearPlane;
    filmPlanDepth = glm::abs(nearPlane / farPlane);
}

void Camera::setFarPlane (float _farPlane) {
    farPlane = _farPlane;
    filmPlanDepth = glm::abs(nearPlane / farPlane);
}

void Camera::setScreenSize (int _screenWidth, int _screenHeight) {
    screenWidth = _screenWidth;
    screenHeight = _screenHeight;
    screenWidthRatio = screenWidth / screenHeight;
}

glm::mat4 Camera::getModelViewMatrix() {
    return glm::transpose(basis) * glm::translate(glm::mat4(1.0), -eyePoint);
}


void Camera::rotateV(float degrees) {
    // y: yaw
    basis = glm::rotate(glm::mat4(1.0), glm::radians(degrees), glm::vec3(basis[1][0], basis[1][1], basis[1][2])) * basis;
}

void Camera::rotateU(float degrees) {
    // x: pitch
    basis = glm::rotate(glm::mat4(1.0), glm::radians(degrees), glm::vec3(basis[0][0], basis[0][1], basis[0][2])) * basis;
}

void Camera::rotateW(float degrees) {
    // z: roll
    basis = glm::rotate(glm::mat4(1.0), glm::radians(degrees), glm::vec3(basis[2][0], basis[2][1], basis[2][2])) * basis;
}

void Camera::translate(glm::vec3 v) {
    eyePoint.x += v.x;
    eyePoint.y += v.y;
    eyePoint.z -= v.z;
}

void Camera::rotate(glm::vec3 point, glm::vec3 axis, float degrees) {
    glm::mat4 I(1.0);
    glm::mat4 tOrigin = glm::translate(I, -point);
    glm::mat4 rotate = glm::rotate(I, glm::radians(degrees), axis);
    glm::mat4 tPoint = glm::translate(I, point);

    eyePoint = tPoint * rotate * tOrigin * glm::vec4(eyePoint.x, eyePoint.y, eyePoint.z, 1.0f);
    basis = rotate * basis;
}

glm::vec3 Camera::getEyePoint() {
	return eyePoint;
}

glm::vec3 Camera::getLookVector() {
	return {-basis[2][0], -basis[2][1], -basis[2][2]};
}

glm::vec3 Camera::getUpVector() {
	return {basis[1][0], basis[1][1], basis[1][2]};
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

float Camera::getFilmPlanDepth() {
    return filmPlanDepth;
}

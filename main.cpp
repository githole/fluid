#include <iostream>
#include <vector>

#include <SDL.h>
#include <Windows.h>
#include <gl\GL.h>
#include <gl\GLU.h>

#include "vec.h"
#include "constant.h"
#include "random.h"


Random rnd(0);

const double SimulationSpaceWidth = 1.0; // m
const double SimulationSpaceHeight = 1.0; // m

struct Particle {
	Vec position_;
	Vec velocity_;
	Vec force_;
	double pressure_;
	double density_;

	std::vector<Particle*> neighbor_;

	Vec PRESSURE_FORCE;
};

void init();
int simulate();

static int fire = 0, del = 0;

int ProcessSDLEvents() {
	SDL_Event eve;
	while (SDL_PollEvent(&eve)) {
		switch(eve.type) {
			
		case SDL_MOUSEBUTTONDOWN:
		if(eve.button.button == SDL_BUTTON_LEFT){
			fire = 1;
		}
		if(eve.button.button == SDL_BUTTON_RIGHT){
			del = 1;
		}
		break;

		case SDL_QUIT:
			return -1;
		}

	}

	return 0;
}

int main(int argc, char* argv[])
{    
	const int width = 960;
	const int height = 960;
	const SDL_VideoInfo* info = NULL;

	if (SDL_Init(SDL_INIT_VIDEO) < 0) {
		return -1;
	}

	info = SDL_GetVideoInfo( );

	if (!info) {
		return -1;
	}
	int bpp = info->vfmt->BitsPerPixel;

	SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 5);
	SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 5);
	SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 5);
	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 16);
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);

	if (SDL_SetVideoMode(width, height, bpp, SDL_OPENGL) == 0) {
		return -1;
	}

	init();

	while (true) {
		if (ProcessSDLEvents() < 0)
			break;


		simulate();
	}



	return 0;
}

void draw_circle(const Vec &p, const double radius) {
	glBegin(GL_LINE_LOOP);
	for (int i = 0; i < 16; i ++) {
		const double a = (i/16.0) * 2.0 * kPI;
		const double x = p.x_ + radius * cos(a);
		const double y = p.y_ + radius * sin(a);
		glVertex2f(x, y);
	}
	glEnd();
}

// パーティクル初期化
const int GridNumX = 32;
const int GridNumY = 32;
const double GridWidth = 0.5; // m
const double GridHeight = 0.5; // m
const double Volume = GridWidth * GridHeight; // m^2	
const double WaterDensity = 100; // kg/(m^2)
const double ParticleMass = WaterDensity * Volume / (GridNumX * GridNumY);
const double EffectiveRadius = 0.02;
const double ParticleRadius = 0.5 * EffectiveRadius;
const double InitDensity = WaterDensity;
const double GasStiffness = 32.0;
const double Viscosity = 4.0;
const Vec Gravity = Vec(0.0, -9.8, 0.0);
std::vector<Particle> particle_buffer[3];
int buffer_index = 0;

// カーネル定数
const double Wpoly6 = 4.0 / (kPI * pow(EffectiveRadius, 8));
const double Gpoly6 = -24.0 / (kPI * pow(EffectiveRadius, 8));
	
const double WSpiky = 10.0 / (kPI * pow(EffectiveRadius, 5));
const double GWSpiky = -30.0 / (kPI * pow(EffectiveRadius, 5));

const double LWVisc = 20.0 / (3.0 * kPI * pow(EffectiveRadius, 5));

void init() {
	for (int i = 0; i < 3; i ++)
		particle_buffer[i].resize(GridNumX * GridNumY);

	for (int ix = 0; ix < GridNumX; ++ix) {
		for (int iy = 0; iy < GridNumY; ++iy) {
			particle_buffer[0][iy * GridNumX + ix].position_ = Vec(0.1 + GridWidth * (double)(ix + 0.5) / GridNumX, 0.1 + GridHeight * (double)(iy + 0.5) / GridNumY);
			particle_buffer[0][iy * GridNumX + ix].velocity_ = Vec(0.0, 0.0);
		}
	}
	particle_buffer[1] = particle_buffer[0];
	particle_buffer[2] = particle_buffer[0];
}

const int NGridX = SimulationSpaceWidth / (EffectiveRadius * 2);
const int NGridY = SimulationSpaceHeight / (EffectiveRadius * 2);

std::vector<Particle*> *grid = new std::vector<Particle*>[NGridX * NGridY];

static double tm = 0.0;
void calc_pos(const Vec &v, int *x, int *y) {
	int ix = v.x_ / SimulationSpaceWidth * NGridX;
	int iy = v.y_ / SimulationSpaceHeight * NGridY;
	if (ix < 0) ix = 0;
	if (ix >= NGridX) ix = NGridX - 1;
	if (iy < 0) iy = 0;
	if (iy >= NGridY) iy = NGridY - 1;

	*x = ix;
	*y = iy;
}

double fract(const double d) {
	return d - int(d);
}

bool check(const int x0, const int y0, const int x1, const int y1) {
	if (x0 == x1)
		return y0 != y1;
	
	return true;
}

int simulate() {

	char title[256];
	sprintf(title, "%f", tm);
	SDL_WM_SetCaption(title, NULL);

	if (fire) {
		fire = 0;
		std::vector<Particle> particle(GridNumX * GridNumY);
		for (int ix = 0; ix < GridNumX; ++ix) {
			for (int iy = 0; iy < GridNumY; ++iy) {
				particle[iy * GridNumX + ix].position_ = Vec(0.1 + GridWidth * (double)(ix + 0.5) / GridNumX, 0.5 + GridHeight * (double)(iy + 0.5) / GridNumY);
				particle[iy * GridNumX + ix].velocity_ = Vec(0.0, 0.0);
			}
		}
		particle_buffer[0].insert(particle_buffer[0].end(), particle.begin(), particle.end()); 
		particle_buffer[1].insert(particle_buffer[1].end(), particle.begin(), particle.end()); 
		particle_buffer[2].insert(particle_buffer[2].end(), particle.begin(), particle.end()); 
	}


	// シミュレーション
	const double dt = 0.002;
//	const double dt = 0.001; // Tait

	tm += dt;

	std::vector<Particle> &new_particle = particle_buffer[buffer_index % 3];
	std::vector<Particle> &two_before_particle = particle_buffer[(buffer_index + 1) % 3];
	std::vector<Particle> &particle = particle_buffer[(buffer_index + 2) % 3];
	
	const int particle_num = new_particle.size();
	if (del) {
		del = 0;

		for (int i = 0; i < particle_num / 2; ++i) {
			particle_buffer[0].pop_back();
			particle_buffer[1].pop_back();
			particle_buffer[2].pop_back();
		}
	}
	// 近傍粒子計算
	for (int i = 0; i < NGridX * NGridY; ++i)
		grid[i].clear();
	for (int i = 0; i < particle_num; ++i) {
		int ix, iy;
		calc_pos(particle[i].position_, &ix, &iy);

		grid[iy * NGridX + ix].push_back(&particle[i]);

		//std::cout << ix << " " << iy << ";";
	}

	for (int i = 0; i < particle_num; ++i) {
		int ix0, iy0;
		calc_pos(particle[i].position_, &ix0, &iy0);

		particle[i].neighbor_.clear();
		for (int oy = -1; oy <= 1; oy ++) {
			for (int ox = -1; ox <= 1; ox ++) {
				const int nx = ox + ix0;
				const int ny = oy + iy0;
				if (nx < 0 || NGridX <= nx)
					continue;
				if (ny < 0 || NGridY <= ny)
					continue;
				
				particle[i].neighbor_.insert(particle[i].neighbor_.end(), grid[ny * NGridX + nx].begin(), grid[ny * NGridX + nx].end());
			}
		}
	}

	// 各粒子位置での密度・圧力値計算
	for (int i = 0; i < particle_num; ++i) {
		particle[i].density_ = 0.0;
		// 密度計算
		const int size = particle[i].neighbor_.size();
		for (int j = 0; j < size; ++j) {
			/*
			const double r = (particle[i].position_ - particle[i].neighbor_[j]->position_).length();
			const double q = EffectiveRadius * EffectiveRadius - r * r;
			if (r >= 0 && r < EffectiveRadius)
				particle[i].density_ += ParticleMass * Wpoly6 * q * q * q;
				*/
			const double r2 = (particle[i].position_ - particle[i].neighbor_[j]->position_).length_squared();
			const double q = EffectiveRadius * EffectiveRadius - r2;
			if (r2 < EffectiveRadius * EffectiveRadius)
				particle[i].density_ += ParticleMass * Wpoly6 * q * q * q;
		}
		// 圧力値計算
		particle[i].pressure_ = GasStiffness * (particle[i].density_ - InitDensity);
		
		// Tait
		/*
		const double Gamma = 7.0;
		const double Cs = 5;
		const double B = InitDensity * Cs * Cs / Gamma;
		particle[i].pressure_ = B * (pow(particle[i].density_ / InitDensity, Gamma)  - 1.0);
		*/
		
	}

	// 圧力項計算
	for (int i = 0; i < particle_num; ++i) {
		const int size = particle[i].neighbor_.size();
//		const double pressure_i = particle[i].pressure_ / pow(particle[i].density_, 2.0);
			const double pressure_i = particle[i].pressure_;
		Vec pressure_force;
		Vec viscosity_force;
		for (int j = 0; j < size; ++j) {
			// 同じ粒子だったら考慮しない
			if (&particle[i] == particle[i].neighbor_[j])
				continue;

			const Vec rij = particle[i].position_ - particle[i].neighbor_[j]->position_;
			const Vec vji = particle[i].neighbor_[j]->velocity_ - particle[i].velocity_;
			const double r = rij.length();
			
			if (r > EffectiveRadius)
				continue;

			const double q = EffectiveRadius - r;
//			const double pressure_j = particle[i].neighbor_[j]->pressure_ / pow(particle[i].neighbor_[j]->density_, 2.0);
				const double pressure_j = particle[i].neighbor_[j]->pressure_;

			Vec kernel;
			kernel = GWSpiky * q * q * rij / r;

			// 圧力項
//			pressure_force = pressure_force + ParticleMass * (pressure_i + pressure_j) * kernel;
			pressure_force = pressure_force + ParticleMass * (pressure_i + pressure_j) / (2.0 * particle[i].neighbor_[j]->density_ ) * kernel;
			particle[i].PRESSURE_FORCE = pressure_force;

			// 粘性
			viscosity_force = viscosity_force + ParticleMass * (vji / particle[i].neighbor_[j]->density_) * LWVisc * q;
		}

		particle[i].force_ = Vec();
//		particle[i].force_  = particle[i].force_  + (-particle[i].density_ * pressure_force);
			particle[i].force_  = particle[i].force_  + (-1.0 * pressure_force);
		particle[i].force_  = particle[i].force_  + Viscosity * viscosity_force;

		// 外力
		particle[i].force_  = particle[i].force_  + Gravity * particle[i].density_;
	}


	// 位置更新
	for (int i = 0; i < particle_num; ++i) {
		new_particle[i].velocity_ = particle[i].velocity_ + dt * particle[i].force_ / particle[i].density_;
		new_particle[i].position_ = particle[i].position_ + dt * new_particle[i].velocity_;

		// シミュレーション空間との衝突判定
		Vec hit;
		Vec repulson;
		const double R = 1000.0;
		const double DAMP = 0.0;
		if (new_particle[i].position_.x_ - ParticleRadius < 0.0) {
			repulson.x_ = R * fabs(new_particle[i].position_.x_ - ParticleRadius) + DAMP * new_particle[i].velocity_.x_;
			new_particle[i].position_.x_ = ParticleRadius;
			hit.x_ = -1.0;
		}
		if (new_particle[i].position_.x_ + ParticleRadius > SimulationSpaceWidth) {
			repulson.x_ = -R * fabs(new_particle[i].position_.x_ + ParticleRadius - SimulationSpaceWidth) - DAMP * new_particle[i].velocity_.x_;
			new_particle[i].position_.x_ = SimulationSpaceWidth - ParticleRadius;
			hit.x_ = 1.0;
		}
			
		if (new_particle[i].position_.y_ - ParticleRadius < 0.0) {
			repulson.y_ = R * fabs(new_particle[i].position_.y_ - ParticleRadius) + DAMP * new_particle[i].velocity_.y_;
			new_particle[i].position_.y_ = ParticleRadius;
			hit.y_ = -1.0;
		}
		if (new_particle[i].position_.y_ + ParticleRadius > SimulationSpaceHeight) {
			repulson.y_ = -R * fabs(new_particle[i].position_.y_ + ParticleRadius - SimulationSpaceHeight) - DAMP * new_particle[i].velocity_.y_;
			new_particle[i].position_.y_ = SimulationSpaceHeight - ParticleRadius;
			hit.y_ = 1.0;
		}

		if (hit.x_ == 0.0 && hit.y_ == 0.0 && hit.z_ == 0.0)
			continue;
		//new_particle[i].velocity_ = new_particle[i].velocity_ + repulson;
		new_particle[i].velocity_ = 0.5 * multiply(hit, new_particle[i].velocity_);
		/*
		normalize(hit);
		new_particle[i].velocity_ = new_particle[i].velocity_ - 0.2 * dot(hit, new_particle[i].velocity_) * hit;
		*/
	}

	buffer_index ++;


	// 描画
	const double scale = 1.5;


	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLineWidth(2);
	for (int i = 0; i < particle_num; i ++) {
		const Vec p = Vec(-scale, -scale) / 2.0 + scale * particle[i].position_;

		double pr = particle[i].pressure_ /1500.0;
		double pg = particle[i].density_ /200.0;

		glColor3f(pr + 0.4, 0.4, 0.4);
		draw_circle(p, scale * ParticleRadius / 2);

		/*
		glColor3f(0.5, 0.0, 0.0);
		glBegin(GL_LINE_LOOP);
		glVertex2f(p.x_, p.y_);
		glVertex2f(p.x_ + particle[i].force_.x_ * 0.0001, p.y_+ particle[i].force_.y_ * 0.0001);
		glEnd();
		
		glColor3f(0.0, 0.5, 0.0);
		glBegin(GL_LINE_LOOP);
		glVertex2f(p.x_, p.y_);
		glVertex2f(p.x_ + particle[i].PRESSURE_FORCE.x_ * 0.0001, p.y_+ particle[i].PRESSURE_FORCE.y_ * 0.0001);
		glEnd();*/
	}
	glColor3f(1.0, 1.0, 1.0);
	
	const double w = scale / 2.0;
	glBegin(GL_LINE_LOOP);
		glVertex2f(-w, w);
		glVertex2f( w, w);
		glVertex2f( w, -w);
		glVertex2f( -w, -w);
	glEnd();
	
	SDL_GL_SwapBuffers();
	return 0;
}
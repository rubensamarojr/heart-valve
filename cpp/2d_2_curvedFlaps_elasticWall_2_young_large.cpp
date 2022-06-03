/**
 * @file 	aortic_valve.cpp
 * @brief 	This is the benchmark test for the 2d valve with 2 curved flaps in contact. T. Wick, Flapping and contact FSI computations with the fluid-solid interface-tracking/interface-capturing technique and mesh adaptivity, Comput. Mech. 53 (2014) 29–43.
 * @author 	Rubens A. Amaro Junior
 */

/**
* @brief SPHinXsys Library.
*/
#include "sphinxsys.h"
/**
* @brief Namespace cite here.
*/
using namespace SPH;

/// RELAXATION AND RELOAD
/// Tag for run particle relaxation for the initial body fitted distribution.
bool particle_relax_ = false;
/// Tag for computation start with relaxed body fitted particles distribution.
bool particle_reload_;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 12.0e-2;									//< Channel length
Real DH = 2.50e-2;									//< Channel height
Real Leaflet_L = 0.750e-2;							//< Leaflet length
Real Leaflet_T = 0.065e-2;						//< Leaflet thickness
Real N = 180;									//< Number of particles for creating leaflet curve
Real Rint = 1.085e-2;									//< leaflet curve radius intern
Real Rext = Rint + Leaflet_T;									//< leaflet curve radius extern
Real resolution_ref = Leaflet_T / 4.0;			//< Global reference particle spacing
Real shiftCorr = - 0.5 * resolution_ref;    //< Correction to avoid fluid-solid overlap at t = 0
Vec2d insert_circle_center_top(3.0e-2, 0.5 * DH);	//< Location of the leaflet curve center
Vec2d insert_circle_center_bottom(3.0e-2, -0.5 * DH + shiftCorr);	//< Location of the leaflet curve center
Real DL_sponge = resolution_ref * 20.0;			//< Sponge region to impose inflow condition
Real BW = 0.2e-2;					//< Boundary width
Real Leaflet_base_length = BW;					//< Length of constrained leaflet
Real Leaflet_dist = 1.85e-2;        //< Leaflet distance from ifnlow

Real DL_extension = DL; // Channel length extension

// Domain bounds of the system.
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -0.5 * DH - BW), Vec2d(DL + BW + DL_extension, 0.5 * DH + BW));



//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0e3;						//< Density.
Real U_f = 1.35e-2;							//< freestream velocity.
Real c_f = 1.0 * 10.0 * U_f;						//< Speed of sound.
Real mu_f = 3.0e-3;							//< Dynamics viscosity.
Real Re = rho0_f * U_f * DH / mu_f;			//< Reynolds number.
//----------------------------------------------------------------------
//	Material properties of the solid.
//----------------------------------------------------------------------
Real rho0_s = 1.2e3;						//< Density.
Real poisson = 0.3;							//< Poisson ratio.
Real Youngs_modulus_left = 1.3e7 * 0.981;				//< Young modulus.
Real Youngs_modulus_right = 1.3e4  * 0.981;				//< Young modulus.

Real gravity_x = 0.0;						//< X-value of gravity.
Real gravity_y = 0.0;						//< Y-value of gravity.

// Initial fluid velocity
Vec3d initial_velocity = Vec3d(0.0, 0.0, 0.0);
//Vec3d initial_velocity = Vec3d(40.0 * 1.0e-4, 0.0, 0.0);
//Vec3d initial_velocity = Vec3d(135.0 * 1.0e-4, 0.0, 0.0);

// Reduce fluid time step
Real dt_reduct = 1.0;

// Final time of simulation
Real Time_Simulation = 1.5;

//----------------------------------------------------------------------
//	Define geometry of SPH bodies
//----------------------------------------------------------------------
// create a water block shape.
std::vector<Vecd> createWaterBlockShape()
{
	//geometry
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(-DL_sponge, -0.5 * DH + shiftCorr));
	water_block_shape.push_back(Vecd(-DL_sponge,  0.5 * DH));
	water_block_shape.push_back(Vecd(DL + DL_extension,  0.5 * DH));
	water_block_shape.push_back(Vecd(DL + DL_extension, -0.5 * DH + shiftCorr));
	water_block_shape.push_back(Vecd(-DL_sponge, -0.5 * DH + shiftCorr));

	return water_block_shape;
}
// create a water block buffer shape.
MultiPolygon createInflowBufferShape()
{
	std::vector<Vecd> inflow_buffer_shape;
	inflow_buffer_shape.push_back(Vecd(-DL_sponge, -0.5  *DH));
	inflow_buffer_shape.push_back(Vecd(-DL_sponge,  0.5 * DH));
	inflow_buffer_shape.push_back(Vecd(0.0,  0.5 * DH));
	inflow_buffer_shape.push_back(Vecd(0.0, -0.5 * DH));
	inflow_buffer_shape.push_back(Vecd(-DL_sponge, -0.5 * DH));

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(inflow_buffer_shape, ShapeBooleanOps::add);
	return multi_polygon;
}

// create outer wall shape.
std::vector<Vecd> createOuterWallShape()
{
	// geometry
	std::vector<Vecd> outer_wall_shape;
	outer_wall_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH - BW));
	outer_wall_shape.push_back(Vecd(-DL_sponge - BW,  0.5 * DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW + DL_extension,  0.5 * DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW + DL_extension, -0.5 * DH - BW));
	outer_wall_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH - BW));

	return outer_wall_shape;
}

// create inner wall shape.
std::vector<Vecd> createInnerWallShape()
{
	std::vector<Vecd> inner_wall_shape;
	inner_wall_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH));
	inner_wall_shape.push_back(Vecd(-DL_sponge - BW,  0.5 * DH));
	inner_wall_shape.push_back(Vecd(DL + BW + DL_extension,  0.5 * DH));
	inner_wall_shape.push_back(Vecd(DL + BW + DL_extension, -0.5 * DH));
	inner_wall_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH));

	return inner_wall_shape;
}

// create top inner wall shape.
std::vector<Vecd> createTopInnerWallShape()
{
	std::vector<Vecd> inner_wall_shape;
	inner_wall_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH - BW));
	inner_wall_shape.push_back(Vecd(-DL_sponge - BW,  0.5 * DH));
	inner_wall_shape.push_back(Vecd(DL + BW + DL_extension,  0.5 * DH));
	inner_wall_shape.push_back(Vecd(DL + BW + DL_extension, -0.5 * DH - BW));
	inner_wall_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH - BW));

	return inner_wall_shape;
}

// create bottom inner wall shape.
std::vector<Vecd> createBottomInnerWallShape()
{
	std::vector<Vecd> inner_wall_shape;
	inner_wall_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH));
	inner_wall_shape.push_back(Vecd(-DL_sponge - BW,  0.5 * DH + BW));
	inner_wall_shape.push_back(Vecd(DL + BW + DL_extension,  0.5 * DH + BW));
	inner_wall_shape.push_back(Vecd(DL + BW + DL_extension, -0.5 * DH));
	inner_wall_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH));

	return inner_wall_shape;
}

// create center wall shape.
std::vector<Vecd> createCentralWallShape()
{
	// geometry
	std::vector<Vecd> outer_wall_shape;
	outer_wall_shape.push_back(Vecd(-DL_sponge, -0.5 * DH - BW));
	outer_wall_shape.push_back(Vecd(-DL_sponge,  0.5 * DH + BW));
	outer_wall_shape.push_back(Vecd(DL,  0.5 * DH + BW));
	outer_wall_shape.push_back(Vecd(DL, -0.5 * DH - BW));
	outer_wall_shape.push_back(Vecd(-DL_sponge, -0.5 * DH - BW));

	return outer_wall_shape;
}


//----------------------------------------------------------------------
// Create top leaflet shape
//----------------------------------------------------------------------
Vec2d TL0(Leaflet_L + insert_circle_center_top[0], insert_circle_center_top[1] - Rint);
Vec2d TL1(Leaflet_L + insert_circle_center_top[0], insert_circle_center_top[1] - Rext);
Vec2d TL2(Leaflet_dist + Leaflet_T, insert_circle_center_top[1]);

std::vector<Vecd> createtopleafletShape()
{
	std::vector<Vecd> leaflet_shape;
	for (int i = 0; i < N + 1; ++i)
	{
			leaflet_shape.push_back(Vecd(insert_circle_center_top[0] - Rint * cos(i * Pi / (2.0 * N)),
			insert_circle_center_top[1] - Rint * sin(i * Pi / (2.0 * N))));
	}
	leaflet_shape.push_back(TL0);
	leaflet_shape.push_back(TL1);
	for (int i = 0; i < N + 1; ++i)
	{
			leaflet_shape.push_back(Vecd(insert_circle_center_top[0] - Rext * sin(i * Pi / (2.0 * N)), 
			insert_circle_center_top[1] - Rext * cos(i * Pi / (2.0 * N))));
	}
	leaflet_shape.push_back(TL2);
	
	return leaflet_shape;
}

//----------------------------------------------------------------------
// Create bottom leaflet shape
//----------------------------------------------------------------------
Vec2d BL0(Leaflet_L + insert_circle_center_bottom[0], insert_circle_center_bottom[1] + Rext);
Vec2d BL1(Leaflet_L + insert_circle_center_bottom[0], insert_circle_center_bottom[1] + Rint);
Vec2d BL2(Leaflet_dist, insert_circle_center_bottom[1]);

std::vector<Vecd> createbottomleafletShape()
{
	std::vector<Vecd> leaflet_shape;
	for (int i = 0; i < N + 1; ++i)
	{
			leaflet_shape.push_back(Vecd(insert_circle_center_bottom[0] - Rext * cos(i * Pi / (2.0 * N)), 
			insert_circle_center_bottom[1] + Rext * sin(i * Pi / (2.0 * N))));
	}
	leaflet_shape.push_back(BL0);
	leaflet_shape.push_back(BL1);
	for (int i = 0; i < N + 1; ++i)
	{
			leaflet_shape.push_back(Vecd(insert_circle_center_bottom[0] - Rint * sin(i * Pi / (2.0 * N)), 
			insert_circle_center_bottom[1] + Rint * cos(i * Pi / (2.0 * N))));
	}
	leaflet_shape.push_back(BL2);
	
	return leaflet_shape;
}

// create the top Leaflet base for constrain.
MultiPolygon createTopLeafletBaseShape()
{
	// Geometry definition
	std::vector<Vecd> outer_wall_shape = createOuterWallShape();
	std::vector<Vecd> inner_wall_shape  = createTopInnerWallShape();
	std::vector<Vecd> water_block_shape = createWaterBlockShape();
	//std::vector<Vecd> top_leaflet_shape = createtopleafletShape();
	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
	multi_polygon.addAPolygon(inner_wall_shape,	 ShapeBooleanOps::sub);
	multi_polygon.addAPolygon(water_block_shape, ShapeBooleanOps::sub);
	//multi_polygon.addAPolygon(top_leaflet_shape,	 ShapeBooleanOps::sub);
	return multi_polygon;
}

// create the bottom Leaflet base for constrain.
MultiPolygon createBottomLeafletBaseShape()
{
	// Geometry definition
	std::vector<Vecd> outer_wall_shape = createOuterWallShape();
	std::vector<Vecd> inner_wall_shape  = createBottomInnerWallShape();
	std::vector<Vecd> water_block_shape = createWaterBlockShape();
	//std::vector<Vecd> bottom_leaflet_shape = createbottomleafletShape();
	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
	multi_polygon.addAPolygon(inner_wall_shape,	 ShapeBooleanOps::sub);
	multi_polygon.addAPolygon(water_block_shape, ShapeBooleanOps::sub);
	//multi_polygon.addAPolygon(bottom_leaflet_shape,	 ShapeBooleanOps::sub);
	return multi_polygon;
}

// create the top Leaflet ends for constrain.
MultiPolygon createTopLeafletEndShape()
{
	// Geometry definition
	std::vector<Vecd> outer_wall_shape = createOuterWallShape();
	std::vector<Vecd> inner_wall_shape  = createTopInnerWallShape();
	std::vector<Vecd> central_wall_shape  = createCentralWallShape();
	std::vector<Vecd> water_block_shape = createWaterBlockShape();
	//std::vector<Vecd> top_leaflet_shape = createtopleafletShape();
	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
	multi_polygon.addAPolygon(inner_wall_shape,	 ShapeBooleanOps::sub);
	multi_polygon.addAPolygon(central_wall_shape, ShapeBooleanOps::sub);
        multi_polygon.addAPolygon(water_block_shape, ShapeBooleanOps::sub);
	//multi_polygon.addAPolygon(top_leaflet_shape,	 ShapeBooleanOps::sub);
	return multi_polygon;
}

// create the bottom Leaflet ends for constrain.
MultiPolygon createBottomLeafletEndShape()
{
	// Geometry definition
	std::vector<Vecd> outer_wall_shape = createOuterWallShape();
	std::vector<Vecd> inner_wall_shape  = createBottomInnerWallShape();
	std::vector<Vecd> central_wall_shape  = createCentralWallShape();
	std::vector<Vecd> water_block_shape = createWaterBlockShape();
	//std::vector<Vecd> top_leaflet_shape = createtopleafletShape();
	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
	multi_polygon.addAPolygon(inner_wall_shape,	 ShapeBooleanOps::sub);
	multi_polygon.addAPolygon(central_wall_shape, ShapeBooleanOps::sub);
        multi_polygon.addAPolygon(water_block_shape, ShapeBooleanOps::sub);
	//multi_polygon.addAPolygon(top_leaflet_shape,	 ShapeBooleanOps::sub);
	return multi_polygon;
}

// Define case dependent bodies material, constraint and boundary conditions.
// Fluid body definition
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &system, std::string body_name)
		: FluidBody(system, body_name)
	{
		// Geometry definition.
		std::vector<Vecd> water_block_shape = createWaterBlockShape();
		std::vector<Vecd> top_leaflet_shape = createtopleafletShape();
		std::vector<Vecd> bottom_leaflet_shape = createbottomleafletShape();
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(water_block_shape, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(top_leaflet_shape, ShapeBooleanOps::sub);
		multi_polygon.addAPolygon(bottom_leaflet_shape, ShapeBooleanOps::sub);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
//----------------------------------------------------------------------
//	Wall boundary body with cases-dependent geometries.
//----------------------------------------------------------------------
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &sph_system, const std::string &body_name)
		: SolidBody(sph_system, body_name)
	{
		// Geometry definition.
		std::vector<Vecd> outer_wall_shape = createOuterWallShape();
		std::vector<Vecd> inner_wall_shape 	= createInnerWallShape();
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(inner_wall_shape,	 ShapeBooleanOps::sub);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};

// Definition of the wall top
class WallBoundaryTop : public SolidBody
{
public:
	WallBoundaryTop(SPHSystem &sph_system, const std::string &body_name)
		: SolidBody(sph_system, body_name)
	{
		// Geometry definition.
		std::vector<Vecd> outer_wall_shape = createOuterWallShape();
		std::vector<Vecd> inner_wall_shape 	= createTopInnerWallShape();
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(inner_wall_shape,	 ShapeBooleanOps::sub);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};

// Definition of the wall bottom
class WallBoundaryBottom : public SolidBody
{
public:
	WallBoundaryBottom(SPHSystem &sph_system, const std::string &body_name)
		: SolidBody(sph_system, body_name)
	{
		// Geometry definition.
		std::vector<Vecd> outer_wall_shape = createOuterWallShape();
		std::vector<Vecd> inner_wall_shape 	= createBottomInnerWallShape();
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(inner_wall_shape,	 ShapeBooleanOps::sub);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};

// Definition of the top inserted body as a elastic structure.
class InsertedBodyTop : public SolidBody
{
public:
	InsertedBodyTop(SPHSystem &system, const std::string body_name)
		: SolidBody(system, body_name, makeShared<SPHAdaptation>(1.15, 2.0))
	{
		// Geometry definition
		std::vector<Vecd> outer_wall_shape = createOuterWallShape();
		std::vector<Vecd> inner_wall_shape 	= createTopInnerWallShape();
		std::vector<Vecd> water_block_shape = createWaterBlockShape();
		std::vector<Vecd> top_leaflet_shape 	= createtopleafletShape();
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(inner_wall_shape,	 ShapeBooleanOps::sub);
		multi_polygon.addAPolygon(water_block_shape, ShapeBooleanOps::sub);
		multi_polygon.addAPolygon(top_leaflet_shape,	 ShapeBooleanOps::add);
		MultiPolygonShape multi_polygon_shape(multi_polygon);
		body_shape_.add<LevelSetShape>(this, multi_polygon_shape);
	}
};

// Definition of the bottom inserted body as a elastic structure.
class InsertedBodyBottom : public SolidBody
{
public:
	InsertedBodyBottom(SPHSystem &system, const std::string body_name)
		: SolidBody(system, body_name, makeShared<SPHAdaptation>(1.15, 2.0))
	{
		// Geometry definition
		std::vector<Vecd> outer_wall_shape = createOuterWallShape();
		std::vector<Vecd> inner_wall_shape 	= createBottomInnerWallShape();
		std::vector<Vecd> water_block_shape = createWaterBlockShape();
		std::vector<Vecd> bottom_leaflet_shape 	= createbottomleafletShape();
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(inner_wall_shape,	 ShapeBooleanOps::sub);
		multi_polygon.addAPolygon(water_block_shape, ShapeBooleanOps::sub);
		multi_polygon.addAPolygon(bottom_leaflet_shape,	 ShapeBooleanOps::add);
		MultiPolygonShape multi_polygon_shape(multi_polygon);
		body_shape_.add<LevelSetShape>(this, multi_polygon_shape);
	}
};

//	Case dependent inflow boundary condition.
class ParabolicInflow : public fluid_dynamics::InflowBoundaryCondition
{
	Real u_ave_, u_ref_, t_ref;

public:
	ParabolicInflow(FluidBody &fluid_body, BodyPartByCell &constrained_region)
		: InflowBoundaryCondition(fluid_body, constrained_region)
		{
			u_ave_ = 0.0;
			u_ref_ = U_f;
			t_ref = 1.0;
		}
		//: InflowBoundaryCondition(fluid_body, constrained_region),
		//u_ave_(0.0), u_ref_(0.13889), t_ref(1.0) {}
		//u_ave_(0.0), u_ref_(0.5), t_ref(1.0) {}

	Vecd getTargetVelocity(Vecd &position, Vecd &velocity) override
	{
		Real u = velocity[0];
		Real v = velocity[1];
		if (position[0] < 0.0)
		{
			//u = (DH / 2.0 + position[1]) * ( DH / 2.0 - position[1]) * u_ave_;
      u = (-6.0 * position[1] * position[1] / (DH * DH) + 1.5) * u_ave_ / 1.5;
      v = 0.0;
		}
		return Vecd(u, v);
	}

	void setupDynamics(Real dt = 0.0) override
	{
		Real run_time = GlobalStaticVariables::physical_time_;
		///u_ave_ = u_ref_ * (sin(2.0 * Pi * run_time ) + 1.1);
		///u_ave_ = 0.5 * (u_ref_ * (sin(2.0 * Pi * run_time + 1.5 * Pi ) + 1.0));
		///u_ave_ *= exp(- 1.0 / 100 * run_time);
		Real time_new = run_time;
		if(run_time > 1.0)
				time_new -= 1.0;
		if(run_time > 2.0)
				time_new -= 2.0;
    
		if(time_new <= 0.20)
		{
			u_ref_ = (-9291.6667 * time_new * time_new + 1893.3333 * time_new + 40.0); 
		}
		else if(time_new <= 0.25)
		{
			u_ref_ = (120.0 * (time_new - 0.2) + 47.0); 
		}
		else if(time_new <= 0.30)
		{
			u_ref_ = (-220.0 * (time_new - 0.25) + 53.0); 
		}
		else if(time_new <= 0.37)
		{
			u_ref_ = (114.2857142 * (time_new - 0.3) + 42.0); 
		}
		else if(time_new <= 0.90)
		{
			u_ref_ = (50.0 * exp(-0.4210255685 * (time_new - 0.37))); 
		}
		else
		{
			u_ref_ = 40.0; 
		}
		
		u_ave_ = u_ref_ * 1.0e-4;
		
	}
};

//application dependent initial condition
class FluidInitialCondition : public fluid_dynamics::FluidInitialCondition
{
public:
	FluidInitialCondition(FluidBody &water) : fluid_dynamics::FluidInitialCondition(water){};
protected:
	void Update(size_t index_i, Real dt) override 
	{
		// initial velocity profile
		vel_n_[index_i][0] = (-6.0 * pos_n_[index_i][1] * pos_n_[index_i][1] / (DH * DH) + 1.5) * initial_velocity[0] / 1.5;
		vel_n_[index_i][1] = initial_velocity[1];
	};
};

class LeafletInitialCondition	: public solid_dynamics::ElasticDynamicsInitialCondition
{
public:
	explicit LeafletInitialCondition(SolidBody &body_name)	: solid_dynamics::ElasticDynamicsInitialCondition(body_name){};

protected:
	void Update(size_t index_i, Real dt) override
	{
		// initial velocity profile
		vel_n_[index_i][0] = (-6.0 * pos_n_[index_i][1] * pos_n_[index_i][1] / (DH * DH) + 1.5) * initial_velocity[0] / 1.5;
		vel_n_[index_i][1] = initial_velocity[1];
	};
};

// Leaflet observer particle generator
class TopLeafletObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	TopLeafletObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		positions_volumes_.push_back(std::make_pair(0.5 * (TL0 + TL1), 0.0));
	}
};

class BottomLeafletObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	BottomLeafletObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		positions_volumes_.push_back(std::make_pair(0.5 * (BL0 + BL1), 0.0));
	}
};
// fluid observer particle generator
class FluidObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	FluidObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		/** A line of measruing points at the entrance of the channel */
		size_t number_observation_points = 21;
		Real range_of_measure = DH - resolution_ref * 4.0;
		Real start_of_measure = resolution_ref * 2.0 - 0.5 * DH;
		/** the measuring particles */
		for (size_t i = 0; i < number_observation_points; ++i) {
			Vec2d point_coordinate(0.0, range_of_measure * Real(i) / Real(number_observation_points - 1) + start_of_measure);
			positions_volumes_.push_back(std::make_pair(point_coordinate, 0.0));
		}
	}
};

int main(int ac, char *av[])
{
	//Tag for computation start with relaxed body fitted particles distribution
	if (particle_relax_ == true)
		particle_reload_ = false;
	else
		particle_reload_ = true;
   
	/**	Build up the environment of a SPHSystem. */
	SPHSystem system(system_domain_bounds, resolution_ref);
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	system.run_particle_relaxation_ = particle_relax_;
	/** Tag for computation start with relaxed body fitted particles distribution. */
	system.reload_particles_ = particle_reload_;
	/** Tag for computation from restart files. 0: start with initial condition. */
	system.restart_step_ = 0;
	/** handle command line arguments */
#ifdef BOOST_AVAILABLE
	system.handleCommandlineOptions(ac, av);
#endif
	/** output environment. */
	In_Output in_output(system);
	//----------------------------------------------------------------------
	/**	Creating body, materials and particles for water block. */
	WaterBlock water_block(system, "WaterBody");
	//WaterMaterial* water_material = new WaterMaterial();
	//FluidParticles fluid_particles(water_block, water_material);
	FluidParticles fluid_particles(water_block, makeShared<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f));

	// Initial fluid velocity field
	FluidInitialCondition fluid_initial_velocity(water_block);

	//WallBoundaryTop wall_boundary_top(system, "WallTop");
	//SolidParticles wall_boundary_top_particles(wall_boundary_top);

	//WallBoundaryBottom wall_boundary_bottom(system, "WallBottom");
	//SolidParticles wall_boundary_bottom_particles(wall_boundary_bottom);

	//----------------------------------------------------------------------
	//	Creating body, materials and particles for inserted body.
	//----------------------------------------------------------------------
	InsertedBodyTop inserted_body_top(system, "InsertedBodyTop");
	InsertedBodyBottom inserted_body_bottom(system, "InsertedBodyBottom");
	
	SharedPtr<ParticleGenerator> inserted_body_top_particle_generator = makeShared<ParticleGeneratorLattice>();
	SharedPtr<ParticleGenerator> inserted_body_bottom_particle_generator = makeShared<ParticleGeneratorLattice>();
	if (!system.run_particle_relaxation_ && system.reload_particles_) {
		inserted_body_top_particle_generator = makeShared<ParticleGeneratorReload>(in_output, inserted_body_top.getBodyName());
		inserted_body_bottom_particle_generator = makeShared<ParticleGeneratorReload>(in_output, inserted_body_bottom.getBodyName());
	}
	
	SharedPtr<LinearElasticSolid> inserted_body_top_material = makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus_right, poisson);
	SharedPtr<LinearElasticSolid> inserted_body_bottom_material = makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus_right, poisson);
	//SharedPtr<NeoHookeanSolid> inserted_body_bottom_material = makeShared<NeoHookeanSolid>(rho0_s, Youngs_modulus_right, poisson);
	//SharedPtr<NeoHookeanSolid> inserted_body_top_material = makeShared<NeoHookeanSolid>(rho0_s, Youngs_modulus_right, poisson);
	ElasticSolidParticles inserted_body_top_particles(inserted_body_top, inserted_body_top_material, inserted_body_top_particle_generator);
	ElasticSolidParticles inserted_body_bottom_particles(inserted_body_bottom, inserted_body_bottom_material, inserted_body_bottom_particle_generator);
	// Initial solid velocity field
	LeafletInitialCondition inserted_body_top_initial_velocity(inserted_body_top);
	LeafletInitialCondition inserted_body_bottom_initial_velocity(inserted_body_bottom);
   
	///WallBoundary wall_boundary(system, "Wall");
	///SolidParticles solid_particles(wall_boundary, makeShared<Solid>(rho0_s, inserted_body_material->ContactStiffness()));

	//----------------------------------------------------------------------
	//	Particle and body creation of leaflet and fluid observers.
	//----------------------------------------------------------------------
	ObserverBody top_leaflet_observer(system, "TopLeafletObserver");
	ObserverParticles top_leaflet_observer_particles(top_leaflet_observer, makeShared<TopLeafletObserverParticleGenerator>());
	ObserverBody bottom_leaflet_observer(system, "BottomLeafletObserver");
	ObserverParticles bottom_leaflet_observer_particles(bottom_leaflet_observer, makeShared<BottomLeafletObserverParticleGenerator>());
	ObserverBody fluid_observer(system, "FluidObserver");
	ObserverParticles flow_observer_particles(fluid_observer, makeShared<FluidObserverParticleGenerator>());
 
	// Add pressure to vtp output
	fluid_particles.addAVariableToWrite<Real>("Pressure");
	//----------------------------------------------------------------------
	//	Define simple data file input and outputs functions.
	//----------------------------------------------------------------------
	//In_Output					in_output(system);
	BodyStatesRecordingToVtp 	write_real_body_states(in_output, system.real_bodies_);
	RestartIO 					restart_io(in_output, system.real_bodies_);
	/** topology */
	BodyRelationInner inserted_body_top_inner(inserted_body_top);
	BodyRelationInner inserted_body_bottom_inner(inserted_body_bottom);
	ComplexBodyRelation water_block_complex(water_block, {&inserted_body_top, &inserted_body_bottom});
	//ComplexBodyRelation water_block_complex(water_block, {&inserted_body_top, &inserted_body_bottom, &wall_boundary_top, &wall_boundary_bottom});
	BodyRelationContact water_body_top_contact(inserted_body_top, {&water_block});
	BodyRelationContact water_body_bottom_contact(inserted_body_bottom, {&water_block});
	BodyRelationContact fluid_observer_contact(fluid_observer, {&water_block});
	BodyRelationContact leaflet_observer_contact_top(top_leaflet_observer, {&inserted_body_top});
	BodyRelationContact leaflet_observer_contact_bottom(bottom_leaflet_observer, {&inserted_body_bottom});
	SolidBodyRelationContact top_bottom_contact(inserted_body_top, {&inserted_body_bottom});
	SolidBodyRelationContact bottom_top_contact(inserted_body_bottom, {&inserted_body_top});
	//SolidBodyRelationContact body_wall_top_contact(inserted_body_top, {&wall_boundary_top});
	//SolidBodyRelationContact body_wall_bottom_contact(inserted_body_bottom, {&wall_boundary_bottom});
	SolidBodyRelationSelfContact inserted_body_top_self_contact(inserted_body_top);
	SolidBodyRelationSelfContact inserted_body_bottom_self_contact(inserted_body_bottom);
	
	//std::cout << "5 " << ".\n";
	/**	Check whether run particle relaxation for body-fitted distribution if chosen. */
	if (system.run_particle_relaxation_)
	{
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the insert body particle position. */
		RandomizePartilePosition random_inserted_body_top_particles(inserted_body_top);
		RandomizePartilePosition random_inserted_body_bottom_particles(inserted_body_bottom);
		/** Write the body state to Vtu file. */
		BodyStatesRecordingToVtp write_inserted_body_top_to_vtp(in_output, {&inserted_body_top});
		BodyStatesRecordingToVtp write_inserted_body_bottom_to_vtp(in_output, {&inserted_body_bottom});
		/** Write the particle reload files. */
		ReloadParticleIO write_particle_reload_files_top(in_output, {&inserted_body_top});
		ReloadParticleIO write_particle_reload_files_bottom(in_output, {&inserted_body_bottom});
   
		/** A Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner_top(inserted_body_top_inner);
		relax_dynamics::RelaxationStepInner relaxation_step_inner_bottom(inserted_body_bottom_inner);
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_inserted_body_top_particles.parallel_exec(0.25);
		random_inserted_body_bottom_particles.parallel_exec(0.25);
		relaxation_step_inner_top.surface_bounding_.parallel_exec();
		relaxation_step_inner_bottom.surface_bounding_.parallel_exec();
		write_inserted_body_top_to_vtp.writeToFile(0);
		write_inserted_body_bottom_to_vtp.writeToFile(0);

		/** relax particles of inserted body */
		int ite_p = 0;
		while (ite_p < 1000)
		{
			relaxation_step_inner_top.parallel_exec();
			relaxation_step_inner_bottom.parallel_exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
				write_inserted_body_top_to_vtp.writeToFile(ite_p);
				write_inserted_body_bottom_to_vtp.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process of the inserted body finish !" << std::endl;
		/** Output results. */
		write_particle_reload_files_top.writeToFile(0);
		write_particle_reload_files_bottom.writeToFile(0);
		return 0;
	}


	/** initial condition */
	//BeamInitialCondition top_initial_velocity(inserted_body_top);
	//BeamInitialCondition bottom_initial_velocity(inserted_body_bottom);
 
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	/** Corrected strong configuration for the elastic insertbody */
	solid_dynamics::CorrectConfiguration inserted_body_top_corrected_configuration_in_strong_form(inserted_body_top_inner);
	solid_dynamics::CorrectConfiguration inserted_body_bottom_corrected_configuration_in_strong_form(inserted_body_bottom_inner);
	/**
	 * @brief Methods used for time stepping.
	 */
	Gravity gravity(Vecd(gravity_x, gravity_y));
	/** Initialize particle acceleration. */
	TimeStepInitialization initialize_a_fluid_step(water_block, gravity);
	/** Evaluation of density by summation approach. */
	fluid_dynamics::DensitySummationComplex update_density_by_summation(water_block_complex);
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize get_fluid_advection_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(water_block);
	/** Pressure relaxation using verlet time stepping. */
	/** Here, we do not use Riemann solver for pressure as the flow is viscous. */
	fluid_dynamics::PressureRelaxationWithWall pressure_relaxation(water_block_complex);
	fluid_dynamics::DensityRelaxationRiemannWithWall density_relaxation(water_block_complex);
	/** Computing viscous acceleration with wall. */
	fluid_dynamics::ViscousAccelerationWithWall viscous_acceleration(water_block_complex);
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityCorrectionComplex transport_velocity_correction(water_block_complex);
	/** viscous acceleration and transport velocity correction can be combined because they are independet dynamics. */
	CombinedInteractionDynamics viscous_acceleration_and_transport_correction(viscous_acceleration, transport_velocity_correction);
	/** Computing vorticity in the flow. */
	//fluid_dynamics::VorticityInner compute_vorticity(water_block_inner);
	fluid_dynamics::VorticityInner compute_vorticity(water_block_complex.inner_relation_);
	/** Inflow boundary condition. */
	//ParabolicInflow parabolic_inflow(water_block, new InflowBuffer(water_block, "Buffer"));
	MultiPolygonShape inflow_buffer_shape(createInflowBufferShape());
	BodyRegionByCell inflow_buffer(water_block, "Buffer", inflow_buffer_shape);
	ParabolicInflow parabolic_inflow(water_block, inflow_buffer);

	/** Periodic BCs in x direction. */
	PeriodicConditionInAxisDirectionUsingCellLinkedList periodic_condition(water_block, 0);
	//----------------------------------------------------------------------
	//	Algorithms of FSI.
	//----------------------------------------------------------------------
	/** Compute the force exerted on solid body due to fluid pressure and viscosity. */
	solid_dynamics::FluidPressureForceOnSolid fluid_pressure_force_on_inserted_body_top(water_body_top_contact);
	solid_dynamics::FluidPressureForceOnSolid fluid_pressure_force_on_inserted_body_bottom(water_body_bottom_contact);
	/** Computing viscous force acting on wall with wall model. */
	solid_dynamics::FluidViscousForceOnSolid fluid_viscous_force_on_inserted_body_top(water_body_top_contact);
	solid_dynamics::FluidViscousForceOnSolid fluid_viscous_force_on_inserted_body_bottom(water_body_bottom_contact);
	/** Compute the average velocity of the inserted body */
	solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration_top(inserted_body_top);
	solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration_bottom(inserted_body_bottom);
	//----------------------------------------------------------------------
	//	Algorithms of solid dynamics.
	//----------------------------------------------------------------------
	/** Compute time step size of elastic solid */
	solid_dynamics::AcousticTimeStepSize inserted_body_top_computing_time_step_size(inserted_body_top);
	solid_dynamics::AcousticTimeStepSize inserted_body_bottom_computing_time_step_size(inserted_body_bottom);
	/** Stress relaxation for the inserted body */
	///solid_dynamics::StressRelaxationFirstHalf	inserted_body_top_stress_relaxation_first_half(inserted_body_top_inner);
	///solid_dynamics::StressRelaxationFirstHalf	inserted_body_bottom_stress_relaxation_first_half(inserted_body_bottom_inner);
	solid_dynamics::KirchhoffStressRelaxationFirstHalf	inserted_body_top_stress_relaxation_first_half(inserted_body_top_inner);
	solid_dynamics::KirchhoffStressRelaxationFirstHalf	inserted_body_bottom_stress_relaxation_first_half(inserted_body_bottom_inner);
	solid_dynamics::StressRelaxationSecondHalf	inserted_body_top_stress_relaxation_second_half(inserted_body_top_inner);
	solid_dynamics::StressRelaxationSecondHalf	inserted_body_bottom_stress_relaxation_second_half(inserted_body_bottom_inner);
	/** Algorithms for solid-solid contact. */
	solid_dynamics::ContactDensitySummation top_bottom_update_contact_density(top_bottom_contact);
	solid_dynamics::ContactForce top_bottom_compute_solid_contact_forces(top_bottom_contact);
	/** Algorithms for solid-solid contact. */
	solid_dynamics::ContactDensitySummation bottom_top_update_contact_density(bottom_top_contact);
	solid_dynamics::ContactForce bottom_top_compute_solid_contact_forces(bottom_top_contact);
  // algorithms for solid self contact
  solid_dynamics::DynamicSelfContactForce inserted_body_top_self_contact_forces(inserted_body_top_self_contact);
  solid_dynamics::DynamicSelfContactForce inserted_body_bottom_self_contact_forces(inserted_body_bottom_self_contact);
	/** Constrain region of the inserted body */
	MultiPolygonShape leaflet_base_shape_top(createTopLeafletEndShape());
	MultiPolygonShape leaflet_base_shape_bottom(createBottomLeafletEndShape());
	BodyRegionByParticle leaflet_base_top(inserted_body_top, "TopLeafletBase", leaflet_base_shape_top);
	BodyRegionByParticle leaflet_base_bottom(inserted_body_bottom, "BottomLeafletBase", leaflet_base_shape_bottom);
	solid_dynamics::ConstrainSolidBodyRegion constrain_leaflet_base_top(inserted_body_top, leaflet_base_top);
	solid_dynamics::ConstrainSolidBodyRegion constrain_leaflet_base_bottom(inserted_body_bottom, leaflet_base_bottom);

	/** Update norm */
	solid_dynamics::UpdateElasticNormalDirection inserted_body_top_update_normal(inserted_body_top);
	solid_dynamics::UpdateElasticNormalDirection inserted_body_bottom_update_normal(inserted_body_bottom);
 
	/** Algorithms for solid-solid contact. */
	///solid_dynamics::ContactDensitySummation inserted_body_wall_update_contact_density(body_wall_contact);
	///solid_dynamics::ContactForce inserted_body_wall_compute_solid_contact_forces(body_wall_contact);
 
	//----------------------------------------------------------------------
	// Write observation data into files.
	//----------------------------------------------------------------------
	RegressionTestTimeAveraged<BodyReducedQuantityRecording<solid_dynamics::TotalViscousForceOnSolid>>
		write_total_viscous_force_on_inserted_body_top(in_output, inserted_body_top);
	RegressionTestTimeAveraged<BodyReducedQuantityRecording<solid_dynamics::TotalViscousForceOnSolid>>
		write_total_viscous_force_on_inserted_body_bottom(in_output, inserted_body_bottom);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		write_top_leaflet_tip_displacement("Position", in_output, leaflet_observer_contact_top);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		write_bottom_leaflet_tip_displacement("Position", in_output, leaflet_observer_contact_bottom);
	ObservedQuantityRecording<Vecd>
		write_fluid_velocity("Velocity", in_output, fluid_observer_contact);

	// Pre-simulation	
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	///bottom_initial_velocity.exec();
	fluid_initial_velocity.exec();
	inserted_body_top_initial_velocity.exec();
	inserted_body_bottom_initial_velocity.exec(); 
	/** initialize cell linked lists for all bodies. */
	system.initializeSystemCellLinkedLists();
	/** periodic condition applied after the mesh cell linked list build up
	  * but before the configuration build up. */
	periodic_condition.update_cell_linked_list_.parallel_exec();
	/** initialize configurations for all bodies. */
	system.initializeSystemConfigurations();
	/** computing surface normal direction for the insert body. */
	inserted_body_top_particles.initializeNormalDirectionFromBodyShape();
	inserted_body_bottom_particles.initializeNormalDirectionFromBodyShape();
	/** computing surface normal direction for walls. */
	//wall_boundary_top_particles.initializeNormalDirectionFromBodyShape();
	//wall_boundary_bottom_particles.initializeNormalDirectionFromBodyShape();
	/** computing linear reproducing configuration for the insert body. */
	inserted_body_top_corrected_configuration_in_strong_form.parallel_exec();
	inserted_body_bottom_corrected_configuration_in_strong_form.parallel_exec();
	///top_initial_velocity.exec();
	///bottom_initial_velocity.exec();
	
	/**
	 * @brief The time stepping starts here. Load restart file if necessary.
	 */
	if (system.restart_step_ != 0) {
		//GlobalStaticVariables::physical_time_ = read_restart_files.readRestartFiles(system.restart_step_);
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.restart_step_);
		
		/** Contact model for top. */
		top_bottom_update_contact_density.parallel_exec();
		top_bottom_compute_solid_contact_forces.parallel_exec();
		/** Contact model for bottom. */
		bottom_top_update_contact_density.parallel_exec();
		bottom_top_compute_solid_contact_forces.parallel_exec();
		
		inserted_body_top_self_contact_forces.parallel_exec();
		inserted_body_bottom_self_contact_forces.parallel_exec();
       
		inserted_body_top.updateCellLinkedList();
		inserted_body_bottom.updateCellLinkedList();
		
		top_bottom_contact.updateConfiguration();
		bottom_top_contact.updateConfiguration();
		
		inserted_body_top_self_contact.updateConfiguration();
		inserted_body_bottom_self_contact.updateConfiguration();
       
		water_block.updateCellLinkedList();
		
		periodic_condition.update_cell_linked_list_.parallel_exec();
		
		/** one need update configuration after periodic condition. */
		water_block_complex.updateConfiguration();
		water_body_top_contact.updateConfiguration();
		water_body_bottom_contact.updateConfiguration();
		inserted_body_top_update_normal.parallel_exec();
		inserted_body_bottom_update_normal.parallel_exec();
	}

	/**	Setup computing and initial conditions. */
	size_t number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = Time_Simulation;			/**< End time. */
	Real D_Time = End_Time / 100.0; /**< time stamps for output. */
	Real Dt = 0.0;					/**< Default advection time step sizes for fluid. */
//	Real Dt = 0.1 * D_Time;
	Real dt = 0.0;					/**< Default acoustic time step sizes for fluid. */
	Real dt_s = 0.0;					/**< Default acoustic time step sizes for solid. */
	Real dt_st = 0.0;					/**< Default acoustic time step sizes for solid top. */
	Real dt_sb = 0.0;					/**< Default acoustic time step sizes for solid bottom. */
	
	size_t inner_ite_dt = 0;
	size_t inner_ite_dt_s = 0;
	/**	Statistics for CPU time */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/** First output before the main loop. */
	write_real_body_states.writeToFile();
	write_top_leaflet_tip_displacement.writeToFile(number_of_iterations);
	write_bottom_leaflet_tip_displacement.writeToFile(number_of_iterations);
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			initialize_a_fluid_step.parallel_exec();
			Dt = dt_reduct * get_fluid_advection_time_step_size.parallel_exec();
//			Dt = SMIN(get_fluid_time_step_size.parallel_exec(), inserted_body_computing_time_step_size.parallel_exec());
			update_density_by_summation.parallel_exec();
			viscous_acceleration_and_transport_correction.parallel_exec(Dt);

			/** FSI for viscous force. */
			fluid_viscous_force_on_inserted_body_top.parallel_exec();
			fluid_viscous_force_on_inserted_body_bottom.parallel_exec();
			/** Update normal direction on elastic body.*/
			inserted_body_top_update_normal.parallel_exec();
			inserted_body_bottom_update_normal.parallel_exec();
			inner_ite_dt = 0;
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				//note that dt needs to sufficiently large to avoid divide zero
				//when computing solid average velocity for FSI
				dt = SMIN(dt_reduct * get_fluid_time_step_size.parallel_exec(), Dt);
//				dt = Dt;
				dt_s = SMIN(inserted_body_top_computing_time_step_size.parallel_exec(), inserted_body_bottom_computing_time_step_size.parallel_exec());
				dt = SMIN(dt_s, dt);
				/** Fluid pressure relaxation */
				pressure_relaxation.parallel_exec(dt);
				/** FSI for pressure force. */
				fluid_pressure_force_on_inserted_body_top.parallel_exec();
				fluid_pressure_force_on_inserted_body_bottom.parallel_exec();
				//fluid_pressure_force_on_inserted_body.parallel_exec();
				/** Fluid density relaxation */
				density_relaxation.parallel_exec(dt);

				/** Solid dynamics. */
				inner_ite_dt_s = 0;
				Real dt_s_sum = 0.0;
				average_velocity_and_acceleration_top.initialize_displacement_.parallel_exec();
				average_velocity_and_acceleration_bottom.initialize_displacement_.parallel_exec();
				
				//inserted_body_wall_update_contact_density.parallel_exec();
				//inserted_body_wall_compute_solid_contact_forces.parallel_exec();
        
				while (dt_s_sum < dt)
				{
					dt_st = SMIN(inserted_body_top_computing_time_step_size.parallel_exec(), dt - dt_s_sum);
					dt_sb = SMIN(inserted_body_bottom_computing_time_step_size.parallel_exec(), dt - dt_s_sum);
					dt_s = SMIN(dt_st, dt_sb);
//					dt_s = inserted_body_top_computing_time_step_size.parallel_exec();
//					dt_s = inserted_body_bottom_computing_time_step_size.parallel_exec();

					///inserted_body_wall_update_contact_density.parallel_exec();
					///inserted_body_wall_compute_solid_contact_forces.parallel_exec();
					///inserted_body_self_contact_forces.parallel_exec();
					///inserted_body.updateCellLinkedList();
					///inserted_body_self_contact.updateConfiguration();
					
					top_bottom_update_contact_density.parallel_exec();
					top_bottom_compute_solid_contact_forces.parallel_exec();
					bottom_top_update_contact_density.parallel_exec();
					bottom_top_compute_solid_contact_forces.parallel_exec();
					
					inserted_body_top_self_contact_forces.parallel_exec();
					inserted_body_bottom_self_contact_forces.parallel_exec();
   				
					inserted_body_top.updateCellLinkedList();
					inserted_body_bottom.updateCellLinkedList();
					
					top_bottom_contact.updateConfiguration();
					bottom_top_contact.updateConfiguration();
					
					inserted_body_top_self_contact.updateConfiguration();
					inserted_body_bottom_self_contact.updateConfiguration();
   
					inserted_body_top_stress_relaxation_first_half.parallel_exec(dt_s);
					inserted_body_bottom_stress_relaxation_first_half.parallel_exec(dt_s);
					constrain_leaflet_base_top.parallel_exec();
					constrain_leaflet_base_bottom.parallel_exec();
					inserted_body_top_stress_relaxation_second_half.parallel_exec(dt_s);
					inserted_body_bottom_stress_relaxation_second_half.parallel_exec(dt_s);
					
					//inserted_body.updateCellLinkedList();
					//body_wall_contact.updateConfiguration();

					dt_s_sum += dt_s;
					inner_ite_dt_s++;
				}
				average_velocity_and_acceleration_top.update_averages_.parallel_exec(dt);
				average_velocity_and_acceleration_bottom.update_averages_.parallel_exec(dt);
        
				//inserted_body_top.updateCellLinkedList();
				//inserted_body_bottom.updateCellLinkedList();
				//body_wall_contact.updateConfiguration();

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				parabolic_inflow.parallel_exec();
				inner_ite_dt++;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "	dt / dt_s = " << inner_ite_dt_s << "\n";
				
				if (number_of_iterations % restart_output_interval == 0 && number_of_iterations != system.restart_step_)
					restart_io.writeToFile(number_of_iterations);
			}
			number_of_iterations++;

			/** Water block configuration and periodic condition. */
			periodic_condition.bounding_.parallel_exec();

			water_block.updateCellLinkedList();
			periodic_condition.update_cell_linked_list_.parallel_exec();
			water_block_complex.updateConfiguration();
			/** one need update configuration after periodic condition. */
			top_bottom_update_contact_density.parallel_exec();
			top_bottom_compute_solid_contact_forces.parallel_exec();
			bottom_top_update_contact_density.parallel_exec();
			bottom_top_compute_solid_contact_forces.parallel_exec();
			
			inserted_body_top_self_contact_forces.parallel_exec();
			inserted_body_bottom_self_contact_forces.parallel_exec();
			
			inserted_body_top.updateCellLinkedList();
			inserted_body_bottom.updateCellLinkedList();
			
			top_bottom_contact.updateConfiguration();
			bottom_top_contact.updateConfiguration();
					
			inserted_body_top_self_contact.updateConfiguration();
			inserted_body_bottom_self_contact.updateConfiguration();
			
			water_body_top_contact.updateConfiguration();
			water_body_bottom_contact.updateConfiguration();
			///body_wall_contact.updateConfiguration();
			
			/** write run-time observation into file */
			write_top_leaflet_tip_displacement.writeToFile(number_of_iterations);
			write_bottom_leaflet_tip_displacement.writeToFile(number_of_iterations);
		}

		tick_count t2 = tick_count::now();
		/** write run-time observation into file */
		compute_vorticity.parallel_exec();
		//if (number_of_iterations % screen_output_interval == 0)
			write_real_body_states.writeToFile();
		write_total_viscous_force_on_inserted_body_top.writeToFile(number_of_iterations);
		fluid_observer_contact.updateConfiguration();
		write_fluid_velocity.writeToFile(number_of_iterations);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}

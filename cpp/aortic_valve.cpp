/**
 * @file 	aortic_valve.cpp
 * @brief 	This is the benchmark test for the aortic valve.
 * @author 	Dong Wu
 */

/**
* @brief SPHinXsys Library.
*/
#include "sphinxsys.h"
/**
* @brief Namespace cite here.
*/
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.16;									/**< Channel length. */
Real DH = 0.02;									/**< Channel height. */
Real R = 0.02;									/**< Sinus cavity radius. */
Real Leaflet_L = 0.024;							/**< Leaflet length. */
Real Leaflet_T = 1.6e-3;						/**< Leaflet thickness. */
Real N = 180;									/**< Particle number for creating sinus cavity. */
Vec2d insert_circle_center(0.5 * DL, 0.5 * DH);	/**< Location of the sinus cavity center. */
Real resolution_ref = Leaflet_T / 4.0;			/**< Global reference particle spacing. */
Real DL_sponge = resolution_ref * 20.0;			/**< Sponge region to impose inflow condition. */
Real BW = resolution_ref * 4.0;					/**< Boundary width. */
Real Leaflet_base_length = BW;					/**< Length of constrained leaflet. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -0.5*DH - BW), Vec2d(DL + BW, 0.5*DH + R + BW));


//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0e3;							//< Density.
Real U_f = 0.13889;								//< freestream velocity.
Real c_f = 10.0 * U_f;							//< Speed of sound.
Real mu_f = 4.3e-3;								//< Dynamics viscosity.
Real Re = rho0_f * U_f * DH / mu_f;				//< Reynolds number.
//----------------------------------------------------------------------
//	Material properties of the solid.
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;							//< Density.
Real poisson = 0.49;							//< Poisson ratio.
Real Youngs_modulus = 1.5e6;					//< Young modulus.

/*
Real rho0_f = 1.0;												//< Density.
Real U_f = 1.0;													//< freestream velocity.
Real c_f = 10.0 * U_f;											//< Speed of sound.
Real Re = 100.0;												//< Reynolds number.
//Real mu_f = rho0_f * U_f * (2.0 * insert_circle_radius) / Re;	//< Dynamics viscosity.
Real mu_f = rho0_f * U_f * DL / Re;								//< Dynamics viscosity.

Real rho0_s = 1.0;	 //< Reference density of gate.
Real poisson = 0.49; //< Poisson ratio.
Real Ae = 7.8e3;	 //< Normalized Youngs Modulus.
Real Youngs_modulus = Ae * rho0_f * U_f * U_f;
*/

Real gravity_g = 0.0;											//< Value of gravity.

//----------------------------------------------------------------------
//	Define geometry of SPH bodies
//----------------------------------------------------------------------
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
	//geometry
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(-DL_sponge, -0.5 * DH));
	water_block_shape.push_back(Vecd(-DL_sponge,  0.5 * DH));
	for (int i = 0; i < N + 1; ++i)
	{
		water_block_shape.push_back(Vecd(insert_circle_center[0] - R * cos(i * Pi / N), 
			insert_circle_center[1] + R * sin(i * Pi / N)));
	}
	water_block_shape.push_back(Vecd(DL,  0.5 * DH));
	water_block_shape.push_back(Vecd(DL, -0.5 * DH));
	water_block_shape.push_back(Vecd(-DL_sponge, -0.5 * DH));

	return water_block_shape;
}
/** create a water block buffer shape. */
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

/** create outer wall shape1. */
std::vector<Vecd> createOuterWallShape1()
{
	// geometry
	std::vector<Vecd> outer_wall_shape1;
	outer_wall_shape1.push_back(Vecd(-DL_sponge - BW, -0.5 * DH - BW));
	outer_wall_shape1.push_back(Vecd(-DL_sponge - BW,  0.5 * DH + BW));
	outer_wall_shape1.push_back(Vecd(DL + BW,  0.5 * DH + BW));
	outer_wall_shape1.push_back(Vecd(DL + BW, -0.5 * DH - BW));
	outer_wall_shape1.push_back(Vecd(-DL_sponge - BW, -0.5 * DH - BW));

	return outer_wall_shape1;
}
/** create outer wall shape2. */
std::vector<Vecd> createOuterWallShape2()
{
	// geometry
	std::vector<Vecd> outer_wall_shape2;
	for (int i = 0; i < N + 1; ++i)
	{
		outer_wall_shape2.push_back(Vecd(insert_circle_center[0] - (R + BW) * cos(i * Pi / N), 
			insert_circle_center[1] + (R + BW) * sin(i * Pi / N)));
	}
	outer_wall_shape2.push_back(Vecd(insert_circle_center[0] - R - BW, insert_circle_center[1]));

	return outer_wall_shape2;
}
/**
 * create inner wall shape.
 */
std::vector<Vecd> createInnerWallShape()
{
	std::vector<Vecd> inner_wall_shape;
	inner_wall_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH));
	inner_wall_shape.push_back(Vecd(-DL_sponge - BW,  0.5 * DH));
	inner_wall_shape.push_back(Vecd(DL + BW,  0.5 * DH));
	inner_wall_shape.push_back(Vecd(DL + BW, -0.5 * DH));
	inner_wall_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH));

	return inner_wall_shape;
}
/**
* @brief create leaflet shape
*/
Vec2d BLB(0.5 * DL - R - (Leaflet_base_length + 0.5 * Leaflet_T) * cos(0.25 * Pi), 
	0.5 * DH + (Leaflet_base_length - 0.5 * Leaflet_T) * sin(0.25 * Pi));
Vec2d BLT(0.5 * DL - R - (Leaflet_base_length - 0.5 * Leaflet_T) * cos(0.25 * Pi), 
	0.5 * DH + (Leaflet_base_length + 0.5 * Leaflet_T) * sin(0.25 * Pi));
Vec2d BRT(0.5 * DL - R + (Leaflet_L + 0.5 * Leaflet_T) * cos(0.25 * Pi), 
	0.5 * DH - (Leaflet_L - 0.5 * Leaflet_T) * sin(0.25 * Pi));
Vec2d BRB(0.5 * DL - R + (Leaflet_L - 0.5 * Leaflet_T) * cos(0.25 * Pi), 
	0.5 * DH - (Leaflet_L + 0.5 * Leaflet_T) * sin(0.25 * Pi));
std::vector<Vecd> createleafletShape()
{
	std::vector<Vecd> leaflet_shape;
	leaflet_shape.push_back(BLB);
	leaflet_shape.push_back(BLT);
	leaflet_shape.push_back(BRT);
	leaflet_shape.push_back(BRB);
	leaflet_shape.push_back(BLB);

	return leaflet_shape;
}
/**	Define case dependent bodies material, constraint and boundary conditions. */
/** Fluid body definition */
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &system, std::string body_name)
		: FluidBody(system, body_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> water_block_shape = createWaterBlockShape();
		std::vector<Vecd> leaflet_shape 	= createleafletShape();
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(water_block_shape, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(leaflet_shape, ShapeBooleanOps::sub);
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
		/** Geometry definition. */
		std::vector<Vecd> outer_wall_shape1 = createOuterWallShape1();
		std::vector<Vecd> outer_wall_shape2 = createOuterWallShape2();
		std::vector<Vecd> inner_wall_shape 	= createInnerWallShape();
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(outer_wall_shape1, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(outer_wall_shape2, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(inner_wall_shape,	 ShapeBooleanOps::sub);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
/** Case-dependent material properties */
/*class WaterMaterial : public WeaklyCompressibleFluid
{
public:
	WaterMaterial() : WeaklyCompressibleFluid(real rho0_f, real c_f, real mu_f)
	{
		rho_0_ 	= rho0_f;
		c_0_ 	= c_f;
		mu_ 	= mu_f;
		// supplementary material parameters derived from basic parameters
		assignDerivedMaterialParameters();
	}
};*/
/** Definition of the inserted body as a elastic structure. */
class InsertedBody : public SolidBody
{
public:
	InsertedBody(SPHSystem &system, const std::string body_name)
		: SolidBody(system, body_name, makeShared<SPHAdaptation>(1.15, 2.0))
	{
		// Geometry definition
		std::vector<Vecd> outer_wall_shape1 = createOuterWallShape1();
		std::vector<Vecd> outer_wall_shape2 = createOuterWallShape2();
		std::vector<Vecd> inner_wall_shape 	= createInnerWallShape();
		std::vector<Vecd> water_block_shape = createWaterBlockShape();
		std::vector<Vecd> leaflet_shape 	= createleafletShape();
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(outer_wall_shape1, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(outer_wall_shape2, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(inner_wall_shape,	 ShapeBooleanOps::sub);
		multi_polygon.addAPolygon(water_block_shape, ShapeBooleanOps::sub);
		multi_polygon.addAPolygon(leaflet_shape,	 ShapeBooleanOps::add);
		//body_shape_.add<MultiPolygonShape>(multi_polygon);
		MultiPolygonShape multi_polygon_shape(multi_polygon);
		body_shape_.add<LevelSetShape>(this, multi_polygon_shape);
	}
};

/** the material for inserted body */
/*class InsertBodyMaterial : public LinearElasticSolid
{
public:
	InsertBodyMaterial() : LinearElasticSolid()
	{
		rho_0_ 	= rho0_s;
		E_0_ 	= Youngs_modulus;
		nu_ 	= poisson;

		assignDerivedMaterialParameters();
	}
};*/

// ERROR 
// BodyPartByParticle(solid_body ... solid_body should be SPHBody
/** constraint the cylinder part of the inserted body. */
/*
class LeafletBase : public BodyPartByParticle
{
public:
	LeafletBase(Solid solid_body, const std::string constrained_region_name)
		: BodyPartByParticle(solid_body, constrained_region_name)
	{
		// Geometry definition
		std::vector<Vecd> outer_wall_shape1 = createOuterWallShape1();
		std::vector<Vecd> outer_wall_shape2 = createOuterWallShape2();
		std::vector<Vecd> inner_wall_shape  = createInnerWallShape();
		std::vector<Vecd> water_block_shape = createWaterBlockShape();
		std::vector<Vecd> leaflet_shape 	= createleafletShape();
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(outer_wall_shape1, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(outer_wall_shape2, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(inner_wall_shape,	 ShapeBooleanOps::sub);
		multi_polygon.addAPolygon(water_block_shape, ShapeBooleanOps::sub);
		multi_polygon.addAPolygon(leaflet_shape,	 ShapeBooleanOps::sub);
		// ERROR body_part_shape_ was not declared
		body_part_particles_.add<MultiPolygonShape>(multi_polygon);
		//body_part_shape_ = new LevelSetComplexShape(this, multi_polygon);

		// Tag the constrained particle
		//tagBodyPart();
		// ERROR no matching function, candidate expects 2 arguments
		//BodyPart();
		TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
		tagParticles(tagging_particle_method);
	}
};
*/
/** create the Leaflet base for constrain . */
MultiPolygon createLeafletBaseShape()
{
	/* Geometry definition */
	std::vector<Vecd> outer_wall_shape1 = createOuterWallShape1();
	std::vector<Vecd> outer_wall_shape2 = createOuterWallShape2();
	std::vector<Vecd> inner_wall_shape  = createInnerWallShape();
	std::vector<Vecd> water_block_shape = createWaterBlockShape();
	std::vector<Vecd> leaflet_shape 	= createleafletShape();
	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(outer_wall_shape1, ShapeBooleanOps::add);
	multi_polygon.addAPolygon(outer_wall_shape2, ShapeBooleanOps::add);
	multi_polygon.addAPolygon(inner_wall_shape,	 ShapeBooleanOps::sub);
	multi_polygon.addAPolygon(water_block_shape, ShapeBooleanOps::sub);
	multi_polygon.addAPolygon(leaflet_shape,	 ShapeBooleanOps::sub);
	return multi_polygon;
}

/** inflow buffer */
/*
class InflowBuffer : public BodyPartByCell
{
public:
	//InflowBuffer(FluidBody* fluid_body, std::string constrained_region_name)
	//	: BodyPartByCell(fluid_body, constrained_region_name)
	InflowBuffer(RealBody fluid_body, const std::string constrained_region_name)
		: BodyPartByCell(fluid_body, constrained_region_name)
	{
		// Geometry definition
		std::vector<Vecd> inflow_buffer_shape = createInflowBufferShape();
		//body_part_shape_ = new ComplexShape(constrained_region_name);
		//body_part_shape_->addAPolygon(inflow_buffer_shape, ShapeBooleanOps::add);
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(inflow_buffer_shape, ShapeBooleanOps::add);
		// ERROR body_part_shape_ was not declared
		cell_linked_list_.add<MultiPolygonShape>(multi_polygon);
		
		// Tag the constrained particle
		//tagBodyPart();
		// ERROR no matching function, candidate expects 2 arguments
		//BodyPart();
		TaggingCellMethod tagging_cell_method = std::bind(&BodyRegionByCell::checkNotFar, this, _1, _2);
		tagCells(tagging_cell_method);
	}
};
*/

/**	Case dependent inflow boundary condition. */

class ParabolicInflow : public fluid_dynamics::InflowBoundaryCondition
{
	Real u_ave_, u_ref_, t_ref;

public:
	ParabolicInflow(FluidBody &fluid_body, BodyPartByCell &constrained_region)
		: InflowBoundaryCondition(fluid_body, constrained_region)
		{
			u_ave_ = 0.0;
			u_ref_ = 0.13889;
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
			///u = (-6.0 * position[1] * position[1] / (DH * DH) + 1.5) * u_ave_;
			//u =  6.0 * u_ave_ * position[1] * (DH - position[1]) / DH / DH;
      
			u = (-6.0 * position[1] * position[1] / (DH * DH) + 1.5) * u_ave_ / 1.5;
			v = 0.0;
		}
		return Vecd(u, v);
	}

	void setupDynamics(Real dt = 0.0) override
	{
		Real run_time = GlobalStaticVariables::physical_time_;
		///u_ave_ = u_ref_ * 0.5 * (1.0 + sin(Pi * run_time / t_ref - 0.5 * Pi));
		//u_ave_ = run_time < t_ref ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref)) : u_ref_;

		Real t1 = 0.85;  // Period
		Real a1 = 0.045; // Amplitude
		Real p1 = 0.0;   // Phase
		Real f1 = a1 * sin(run_time * Pi / t1 + p1);

		Real l2 = 1.0;   // Lambda
		Real a2 = 12.0;  // Amplitude
		Real p2 = 0.0;   // Phase
		Real f2 = a2 * (exp(run_time * l2 + p2) - 1.0);

		Real l3 = 0.9;   // Lambda
		Real a3 = -13.0; // Amplitude
		Real p3 = 0.24;  // Phase
		Real f3 = a3 * log10(run_time * l3 + p3);

		Real f4 = a1 * sin((run_time - 1.0) * Pi / t1 + p1);

		if(run_time < 0.28)
		{
			u_ave_ = f1 * f2;
		}
		else if(run_time < 0.85)
		{
			u_ave_ = f1 * f3;
		}
		else if(run_time < 1.0)
		{
			u_ave_ = 0.0;
		}
		else if(run_time < 1.28)
		{
			Real f5 = a2 * (exp((run_time - 1.0) * l2 + p2) - 1.0);
			u_ave_ = f4 * f5;
		}
		else if(run_time < 1.85)
		{
			Real f6 = a3 * log10((run_time - 1.0) * l3 + p3);
			u_ave_ = f4 * f6;
		}
		else
		{
			u_ave_ = 0.0;
		}
    
	}
};


/**	Observer body. */
/*
class LeafletObserver : public FictitiousBody
{
public:
	LeafletObserver(SPHSystem& system, std::string body_name)
		: FictitiousBody(system, body_name, makeShared<SPHAdaptation>(1.15, 1))
	{
		// the measuring particle with zero volume
		// error: ‘body_input_points_volumes_’ was not declared in this scope
		//body_input_points_volumes_.push_back(std::make_pair(0.5 * (BRT + BRB), 0.0));
		positions_volumes_.push_back(std::make_pair(0.5 * (BRT + BRB), 0.0));
	}
};

//	An observer body to measure the flow profile
class FluidObserver : public FictitiousBody
{
public:
	FluidObserver(SPHSystem& system, std::string body_name)
		: FictitiousBody(system, body_name, makeShared<SPHAdaptation>(1.15, 1))
	{
		// A line of measruing points at the entrance of the channel
		size_t number_observation_points = 21;
		Real range_of_measure = DH - resolution_ref * 4.0;
		Real start_of_measure = resolution_ref * 2.0 - 0.5 * DH;
		// the measuring particles
		for (size_t i = 0; i < number_observation_points; ++i) {
			Vec2d point_coordinate(0.0, range_of_measure * Real(i) / Real(number_observation_points - 1) + start_of_measure);
			// error: ‘body_input_points_volumes_’ was not declared in this scope
			//body_input_points_volumes_.push_back(std::make_pair(point_coordinate, 0.0));
			positions_volumes_.push_back(std::make_pair(point_coordinate, 0.0));
		}
	}
};
*/
class LeafletObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	LeafletObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		positions_volumes_.push_back(std::make_pair(0.5 * (BRT + BRB), 0.0));
	}
};
/** fluid observer particle generator */
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
	/**	Build up the environment of a SPHSystem. */
	SPHSystem system(system_domain_bounds, resolution_ref);
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	system.run_particle_relaxation_ = false;
	/** Tag for computation start with relaxed body fitted particles distribution. */
	system.reload_particles_ = true;
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

	//WallBoundary wall_boundary(system, "Wall");
	//SolidParticles wall_particles(wall_boundary);

	//----------------------------------------------------------------------
	//	Creating body, materials and particles for inserted body.
	//----------------------------------------------------------------------
	//InsertedBody* inserted_body = new InsertedBody(system, "InsertedBody");
	InsertedBody inserted_body(system, "InsertedBody");
	//InsertBodyMaterial* insert_body_material = new InsertBodyMaterial();
	//ElasticSolidParticles inserted_body_particles(inserted_body, insert_body_material);
	SharedPtr<ParticleGenerator> inserted_body_particle_generator = makeShared<ParticleGeneratorLattice>();
	if (!system.run_particle_relaxation_ && system.reload_particles_) {
		inserted_body_particle_generator = makeShared<ParticleGeneratorReload>(in_output, inserted_body.getBodyName());
	}
	ElasticSolidParticles inserted_body_particles(inserted_body, 
		makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson), inserted_body_particle_generator);
	//----------------------------------------------------------------------
	//	Particle and body creation of leaflet and fluid observers.
	//----------------------------------------------------------------------
	//LeafletObserver* leaflet_observer = new LeafletObserver(system, "LeafletObserver");
	//BaseParticles leaflet_observer_particles(leaflet_observer);

	ObserverBody leaflet_observer(system, "LeafletObserver");
	ObserverParticles leaflet_observer_particles(leaflet_observer, makeShared<LeafletObserverParticleGenerator>());
	//FluidObserver* fluid_observer = new FluidObserver(system, "FluidObserver");
	ObserverBody fluid_observer(system, "FluidObserver");
	ObserverParticles flow_observer_particles(fluid_observer, makeShared<FluidObserverParticleGenerator>());
	//----------------------------------------------------------------------
	//	Define simple data file input and outputs functions.
	//----------------------------------------------------------------------
	//In_Output					in_output(system);
	BodyStatesRecordingToVtp 	write_real_body_states(in_output, system.real_bodies_);
	//WriteRestart				write_restart_files(in_output, system.real_bodies_);
	//ReadRestart					read_restart_files(in_output, system.real_bodies_);
	RestartIO 					restart_io(in_output, system.real_bodies_);
	/** topology */
	//InnerBodyRelation* water_block_inner			= new InsertBodyRelation(water_block);
	//InnerBodyRelation* inserted_body_inner			= new InsertBodyRelation(inserted_body);
	//ComplexBodyRelation* water_block_complex		= new ComplexBodyRelation(water_block_inner, {&inserted_body});
	//ContactBodyRelation* inserted_body_contact		= new ContactBodyRelation(inserted_body, {&water_block});
	//ContactBodyRelation* leaflet_observer_contact	= new ContactBodyRelation(leaflet_observer, {&inserted_body});
	//ContactBodyRelation* fluid_observer_contact		= new ContactBodyRelation(fluid_observer, {&water_block});

	BodyRelationInner inserted_body_inner(inserted_body);
	//ComplexBodyRelation water_block_complex(water_block, {&wall_boundary, &inserted_body});
	ComplexBodyRelation water_block_complex(water_block, {&inserted_body});
	BodyRelationContact inserted_body_contact(inserted_body, {&water_block});
	BodyRelationContact leaflet_observer_contact(leaflet_observer, {&inserted_body});
	BodyRelationContact fluid_observer_contact(fluid_observer, {&water_block});

	//std::cout << "5 " << ".\n";
	/**	Check whether run particle relaxation for body-fitted distribution if chosen. */
	if (system.run_particle_relaxation_)
	{
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the insert body particle position. */
		RandomizePartilePosition random_inserted_body_particles(inserted_body);
		/** Write the body state to Vtu file. */
		BodyStatesRecordingToVtp write_inserted_body_to_vtp(in_output, {&inserted_body});
		/** Write the particle reload files. */
		ReloadParticleIO write_particle_reload_files(in_output, {&inserted_body});

		/** A Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner(inserted_body_inner);
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_inserted_body_particles.parallel_exec(0.25);
		relaxation_step_inner.surface_bounding_.parallel_exec();
		write_inserted_body_to_vtp.writeToFile(0);

		/** relax particles of inserted body */
		int ite_p = 0;
		while (ite_p < 1000)
		{
			relaxation_step_inner.parallel_exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
				write_inserted_body_to_vtp.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process of the inserted body finish !" << std::endl;
		/** Output results. */
		write_particle_reload_files.writeToFile(0);
		return 0;
	}

	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	/** Corrected strong configuration for the elastic insertbody */
	solid_dynamics::CorrectConfiguration inserted_body_corrected_configuration_in_strong_form(inserted_body_inner);
	/**
	 * @brief Methods used for time stepping.
	 */
	Gravity gravity(Vecd(0.0, -gravity_g));
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
	solid_dynamics::FluidPressureForceOnSolid fluid_pressure_force_on_inserted_body(inserted_body_contact);
	/** Computing viscous force acting on wall with wall model. */
	solid_dynamics::FluidViscousForceOnSolid fluid_viscous_force_on_inserted_body(inserted_body_contact);
	/** Compute the average velocity of the inserted body */
	solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(inserted_body);
	//----------------------------------------------------------------------
	//	Algorithms of solid dynamics.
	//----------------------------------------------------------------------
	/** Compute time step size of elastic solid */
	solid_dynamics::AcousticTimeStepSize inserted_body_computing_time_step_size(inserted_body);
	/** Stress relaxation for the inserted body */
	solid_dynamics::StressRelaxationFirstHalf	inserted_body_stress_relaxation_first_half(inserted_body_inner);
	solid_dynamics::StressRelaxationSecondHalf	inserted_body_stress_relaxation_second_half(inserted_body_inner);
	/** Constrain region of the inserted body */
	//solid_dynamics::ConstrainSolidBodyRegion constrain_leaflet_base(inserted_body, new LeafletBase(inserted_body, "LeafletBase"));
	MultiPolygonShape Leaflet_base_shape(createLeafletBaseShape());
	BodyRegionByParticle leaflet_base(inserted_body, "LeafletBase", Leaflet_base_shape);
	solid_dynamics::ConstrainSolidBodyRegion constrain_leaflet_base(inserted_body, leaflet_base);

	/** Update norm */
	solid_dynamics::UpdateElasticNormalDirection inserted_body_update_normal(inserted_body);
	//----------------------------------------------------------------------
	// Write observation data into files.
	//----------------------------------------------------------------------
	//WriteTotalViscousForceOnSolid write_total_viscous_force_on_inserted_body(in_output, inserted_body);
	//WriteAnObserverQuantity<indexVector, Vecd> write_leaflet_tip_displacement("Position", in_output, leaflet_observer_contact);
	//WriteAnObserverQuantity<indexVector, Vecd> write_fluid_velocity("Velocity", in_output, fluid_observer_contact);

	RegressionTestTimeAveraged<BodyReducedQuantityRecording<solid_dynamics::TotalViscousForceOnSolid>>
		write_total_viscous_force_on_inserted_body(in_output, inserted_body);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		write_leaflet_tip_displacement("Position", in_output, leaflet_observer_contact);
	ObservedQuantityRecording<Vecd>
		write_fluid_velocity("Velocity", in_output, fluid_observer_contact);

	/**
	 * @brief Pre-simulation
	 */
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	/** Using relaxed particle distribution if needed. */
	// error: ‘ReadReloadParticle’ was not declared in this scope
	/*if (system.reload_particles_)
	{
		std::unique_ptr<ReadReloadParticle> reload_inserted_body_particles(new ReadReloadParticle(in_output, {&inserted_body}, {"InsertedBody"}));
		reload_inserted_body_particles.ReadFromFile();
	}*/

	/** initialize cell linked lists for all bodies. */
	system.initializeSystemCellLinkedLists();
	/** periodic condition applied after the mesh cell linked list build up
	  * but before the configuration build up. */
	periodic_condition.update_cell_linked_list_.parallel_exec();
	/** initialize configurations for all bodies. */
	system.initializeSystemConfigurations();
	/** computing surface normal direction for the insert body. */
	inserted_body_particles.initializeNormalDirectionFromBodyShape();
	/** computing linear reproducing configuration for the insert body. */
	inserted_body_corrected_configuration_in_strong_form.parallel_exec();
	
	/**
	 * @brief The time stepping starts here. Load restart file if necessary.
	 */
	if (system.restart_step_ != 0) {
		//GlobalStaticVariables::physical_time_ = read_restart_files.readRestartFiles(system.restart_step_);
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.restart_step_);
		inserted_body.updateCellLinkedList();
		water_block.updateCellLinkedList();
		periodic_condition.update_cell_linked_list_.parallel_exec();
		/** one need update configuration after periodic condition. */
		water_block_complex.updateConfiguration();
		inserted_body_contact.updateConfiguration();
		inserted_body_update_normal.parallel_exec();
	}

	/**	Setup computing and initial conditions. */
	size_t number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 2.0;			/**< End time. */
	Real D_Time = End_Time / 100.0; /**< time stamps for output. */
	Real Dt = 0.0;					/**< Default advection time step sizes for fluid. */
	Real dt = 0.0;					/**< Default acoustic time step sizes for fluid. */
	Real dt_s = 0.0;					/**< Default acoustic time step sizes for solid. */
	size_t inner_ite_dt = 0;
	size_t inner_ite_dt_s = 0;
	/**	Statistics for CPU time */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/** First output before the main loop. */
	write_real_body_states.writeToFile();
	write_leaflet_tip_displacement.writeToFile(number_of_iterations);
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
			Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_density_by_summation.parallel_exec();
			viscous_acceleration_and_transport_correction.parallel_exec(Dt);

			/** FSI for viscous force. */
			fluid_viscous_force_on_inserted_body.parallel_exec();
			/** Update normal direction on elastic body.*/
			inserted_body_update_normal.parallel_exec();
			inner_ite_dt = 0;
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				dt = SMIN(get_fluid_time_step_size.parallel_exec(), Dt);
				/** Fluid pressure relaxation */
				pressure_relaxation.parallel_exec(dt);
				/** FSI for pressure force. */
				fluid_pressure_force_on_inserted_body.parallel_exec();
				/** Fluid density relaxation */
				density_relaxation.parallel_exec(dt);

				/** Solid dynamics. */
				inner_ite_dt_s = 0;
				Real dt_s_sum = 0.0;
				average_velocity_and_acceleration.initialize_displacement_.parallel_exec();
				while (dt_s_sum < dt)
				{
					dt_s = SMIN(inserted_body_computing_time_step_size.parallel_exec(), dt - dt_s_sum);
					inserted_body_stress_relaxation_first_half.parallel_exec(dt_s);
					constrain_leaflet_base.parallel_exec();
					inserted_body_stress_relaxation_second_half.parallel_exec(dt_s);
					dt_s_sum += dt_s;
					inner_ite_dt_s++;
				}
				average_velocity_and_acceleration.update_averages_.parallel_exec(dt);

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
			inserted_body.updateCellLinkedList();
			inserted_body_contact.updateConfiguration();
			/** write run-time observation into file */
			write_leaflet_tip_displacement.writeToFile(number_of_iterations);
		}

		tick_count t2 = tick_count::now();
		/** write run-time observation into file */
		compute_vorticity.parallel_exec();
		write_real_body_states.writeToFile();
		write_total_viscous_force_on_inserted_body.writeToFile(number_of_iterations);
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

// Constants
const G = 39.47841760435743; // Gravitational constant (AU³/M☉·yr²)
const c = 63239.7263;   // Speed of light (AU/year)
const dt = 0.001;       // Time step (years)
const SCALE = 200;      // Visualization scale (pixels per AU)

// Vector3D class
class Vector3D {
    constructor(x, y, z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    add(v) {
        return new Vector3D(this.x + v.x, this.y + v.y, this.z + v.z);
    }

    subtract(v) {
        return new Vector3D(this.x - v.x, this.y - v.y, this.z - v.z);
    }

    multiply(scalar) {
        return new Vector3D(this.x * scalar, this.y * scalar, this.z * scalar);
    }

    dot(v) {
        return this.x * v.x + this.y * v.y + this.z * v.z;
    }

    cross(v) {
        return new Vector3D(
            this.y * v.z - this.z * v.y,
            this.z * v.x - this.x * v.z,
            this.x * v.y - this.y * v.x
        );
    }

    magnitude() {
        return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
    }

    normalize() {
        const m = this.magnitude();
        return m === 0 ? new Vector3D(0, 0, 0) : new Vector3D(this.x / m, this.y / m, this.z / m);
    }
}

// Canvas setup
const canvas = document.getElementById('solarSystem');
const ctx = canvas.getContext('2d');
canvas.width = window.innerWidth * 0.8;
canvas.height = window.innerHeight * 0.8;

// Celestial body class
class CelestialBody {
    constructor(name, mass, radius, position, velocity, color, isSun = false) {
        this.name = name;
        this.mass = mass;
        this.radius = radius;
        this.position = position;
        this.velocity = velocity;
        this.color = color;
        this.isSun = isSun;
        this.trail = [];
        this.maxTrailLength = 1000;
    }

    draw() {
        // Draw trail
        ctx.beginPath();
        ctx.strokeStyle = this.color;
        ctx.lineWidth = 1;
        for (let i = 0; i < this.trail.length - 1; i++) {
            const screenX1 = this.trail[i].x * SCALE + canvas.width / 2;
            const screenY1 = this.trail[i].z * SCALE + canvas.height / 2;
            const screenX2 = this.trail[i + 1].x * SCALE + canvas.width / 2;
            const screenY2 = this.trail[i + 1].z * SCALE + canvas.height / 2;
            ctx.moveTo(screenX1, screenY1);
            ctx.lineTo(screenX2, screenY2);
        }
        ctx.stroke();

        // Draw body
        const screenX = this.position.x * SCALE + canvas.width / 2;
        const screenY = this.position.z * SCALE + canvas.height / 2;
        ctx.beginPath();
        ctx.fillStyle = this.color;
        ctx.arc(
            screenX,
            screenY,
            Math.max(3, this.radius * SCALE),
            0,
            2 * Math.PI
        );
        ctx.fill();

        // Draw name
        ctx.fillStyle = 'white';
        ctx.font = '12px Arial';
        ctx.fillText(this.name, screenX + this.radius * SCALE + 5, screenY);
    }

    updateTrail() {
        this.trail.push(new Vector3D(this.position.x, this.position.y, this.position.z));
        if (this.trail.length > this.maxTrailLength) {
            this.trail.shift();
        }
    }
}

// Light ray class
class LightRay {
    constructor(position, velocity) {
        this.position = position;
        this.velocity = velocity.normalize().multiply(c);
        this.path = [new Vector3D(position.x, position.y, position.z)];
    }

    reset(startX, z) {
        this.position = new Vector3D(startX, 0, z);
        this.velocity = new Vector3D(c, 0, 0);
        this.path = [new Vector3D(this.position.x, this.position.y, this.position.z)];
    }
}

// Initialize solar system
const bodies = [
    // Sun
    new CelestialBody(
        "Sun",
        1.0,
        0.04,
        new Vector3D(0, 0, 0),
        new Vector3D(0, 0, 0),
        "#FFD700",
        true
    ),
    // Mercury
    new CelestialBody(
        "Mercury",
        1.651e-7,
        0.015,
        new Vector3D(0.387, 0, 0),
        new Vector3D(0, 10.06, 0),
        "#A0522D"
    ),
    // Venus
    new CelestialBody(
        "Venus",
        2.447e-6,
        0.02,
        new Vector3D(0.723, 0, 0),
        new Vector3D(0, 7.39, 0),
        "#DEB887"
    ),
    // Earth
    new CelestialBody(
        "Earth",
        3.003e-6,
        0.02,
        new Vector3D(1, 0, 0),
        new Vector3D(0, 6.283, 0),
        "#4169E1"
    ),
    // Mars
    new CelestialBody(
        "Mars",
        3.213e-7,
        0.018,
        new Vector3D(1.524, 0, 0),
        new Vector3D(0, 5.09, 0),
        "#CD5C5C"
    )
];

// Add moons and asteroids
function addMoonsAndAsteroids() {
    // Earth's Moon
    const earth = bodies.find(b => b.name === "Earth");
    if (earth) {
        const moonDist = 0.05;
        const moonSpeed = Math.sqrt(G * earth.mass / moonDist);
        bodies.push(new CelestialBody(
            "Moon",
            3.7e-8,
            0.01,
            earth.position.add(new Vector3D(moonDist, 0, 0)),
            earth.velocity.add(new Vector3D(0, moonSpeed, 0)),
            "#FFFFFF"
        ));
    }

    // Asteroid belt
    const sun = bodies.find(b => b.name === "Sun");
    if (sun) {
        for (let i = 0; i < 100; i++) {
            const angle = Math.random() * 2 * Math.PI;
            const r = 2.2 + Math.random() * (3.2 - 2.2);
            const x = r * Math.cos(angle);
            const z = r * Math.sin(angle);
            const orbitalSpeed = Math.sqrt(G * sun.mass / r);
            const vx = -orbitalSpeed * Math.sin(angle);
            const vz = orbitalSpeed * Math.cos(angle);
            
            bodies.push(new CelestialBody(
                `Asteroid${i}`,
                1e-9,
                0.005,
                new Vector3D(x, 0, z),
                new Vector3D(vx, 0, vz),
                "#808080"
            ));
        }
    }
}

// Initialize light rays
const lightRays = [];
function initLightRays() {
    const numRays = 20;
    const startX = -5;
    const minZ = -3, maxZ = 3;
    
    for (let i = 0; i < numRays; i++) {
        const z = minZ + i * (maxZ - minZ) / (numRays - 1);
        lightRays.push(new LightRay(
            new Vector3D(startX, 0, z),
            new Vector3D(c, 0, 0)
        ));
    }
}

// Calculate acceleration with GR corrections
function calculateAcceleration(body, position) {
    let acc = new Vector3D(0, 0, 0);
    
    for (const other of bodies) {
        if (other === body) continue;
        
        const r = position.subtract(other.position);
        const dist = r.magnitude();
        if (dist < 1e-6) continue;
        
        let factor = G * other.mass / (dist * dist * dist);
        
        // GR correction for Sun's gravity
        if (other.isSun && !body.isSun) {
            const relVel = body.velocity.subtract(other.velocity);
            const l = r.cross(relVel).magnitude();
            const correction = 1 + 3 * (l * l) / (dist * dist * c * c);
            factor *= correction;
        }
        
        acc = acc.subtract(r.multiply(factor));
    }
    
    return acc;
}

// RK4 integration step
function rk4Step(body, dt) {
    const p0 = body.position;
    const v0 = body.velocity;
    
    const a1 = calculateAcceleration(body, p0);
    const k1v = a1.multiply(dt);
    const k1p = v0.multiply(dt);
    
    const v2 = v0.add(k1v.multiply(0.5));
    const p2 = p0.add(k1p.multiply(0.5));
    const a2 = calculateAcceleration(body, p2);
    const k2v = a2.multiply(dt);
    const k2p = v2.multiply(dt);
    
    const v3 = v0.add(k2v.multiply(0.5));
    const p3 = p0.add(k2p.multiply(0.5));
    const a3 = calculateAcceleration(body, p3);
    const k3v = a3.multiply(dt);
    const k3p = v3.multiply(dt);
    
    const v4 = v0.add(k3v);
    const p4 = p0.add(k3p);
    const a4 = calculateAcceleration(body, p4);
    const k4v = a4.multiply(dt);
    const k4p = v4.multiply(dt);
    
    const dp = k1p.add(k2p.multiply(2)).add(k3p.multiply(2)).add(k4p).multiply(1/6);
    const dv = k1v.add(k2v.multiply(2)).add(k3v.multiply(2)).add(k4v).multiply(1/6);
    
    body.position = body.position.add(dp);
    body.velocity = body.velocity.add(dv);
}

// Draw space-time grid
function drawGrid() {
    ctx.strokeStyle = 'rgba(50, 50, 50, 0.5)';
    ctx.lineWidth = 1;
    
    const gridSpacing = 0.5;
    const simWidth = canvas.width / SCALE;
    const simHeight = canvas.height / SCALE;
    
    // Vertical lines
    for (let x = -simWidth/2; x <= simWidth/2; x += gridSpacing) {
        ctx.beginPath();
        for (let y = -simHeight/2; y <= simHeight/2; y += 0.1) {
            const pos = new Vector3D(x, 0, y);
            const bent = bendSpaceTime(pos);
            const screenX = bent.x * SCALE + canvas.width/2;
            const screenY = bent.z * SCALE + canvas.height/2;
            
            if (y === -simHeight/2) {
                ctx.moveTo(screenX, screenY);
            } else {
                ctx.lineTo(screenX, screenY);
            }
        }
        ctx.stroke();
    }
    
    // Horizontal lines
    for (let y = -simHeight/2; y <= simHeight/2; y += gridSpacing) {
        ctx.beginPath();
        for (let x = -simWidth/2; x <= simWidth/2; x += 0.1) {
            const pos = new Vector3D(x, 0, y);
            const bent = bendSpaceTime(pos);
            const screenX = bent.x * SCALE + canvas.width/2;
            const screenY = bent.z * SCALE + canvas.height/2;
            
            if (x === -simWidth/2) {
                ctx.moveTo(screenX, screenY);
            } else {
                ctx.lineTo(screenX, screenY);
            }
        }
        ctx.stroke();
    }
}

// Calculate space-time bending
function bendSpaceTime(point) {
    let offset = new Vector3D(0, 0, 0);
    const k = 0.1;
    
    for (const body of bodies) {
        const d = point.subtract(body.position);
        let r = d.magnitude();
        const minR = 0.2;
        if (r < minR) r = minR;
        
        offset = offset.add(d.normalize().multiply(k * body.mass / r));
    }
    
    return point.subtract(offset);
}

// Animation state
let isRunning = false;
let animationFrameId = null;

// Initialize simulation
addMoonsAndAsteroids();
initLightRays();

// Animation loop
function animate() {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    
    // Draw space-time grid
    drawGrid();
    
    // Update and draw light rays
    for (const ray of lightRays) {
        // Update ray position
        const acc = calculateAcceleration(ray, ray.position);
        ray.velocity = ray.velocity.add(acc.multiply(dt));
        ray.position = ray.position.add(ray.velocity.multiply(dt));
        ray.path.push(new Vector3D(ray.position.x, ray.position.y, ray.position.z));
        
        // Draw ray path
        ctx.beginPath();
        ctx.strokeStyle = 'rgba(0, 255, 255, 0.3)';
        for (let i = 0; i < ray.path.length - 1; i++) {
            const p1 = ray.path[i];
            const p2 = ray.path[i + 1];
            ctx.moveTo(p1.x * SCALE + canvas.width/2, p1.z * SCALE + canvas.height/2);
            ctx.lineTo(p2.x * SCALE + canvas.width/2, p2.z * SCALE + canvas.height/2);
        }
        ctx.stroke();
        
        // Reset ray if out of bounds
        if (ray.position.x > 8) {
            ray.reset(-5, ray.position.z);
        }
        
        // Limit path length
        if (ray.path.length > 500) {
            ray.path.shift();
        }
    }

    // Update positions using RK4
    for (const body of bodies) {
        rk4Step(body, dt);
        body.updateTrail();
        body.draw();
    }

    // Update simulation stats
    updateSimulationStats();

    if (isRunning) {
        animationFrameId = requestAnimationFrame(animate);
    }
}

// Update simulation statistics
function updateSimulationStats() {
    const earth = bodies.find(body => body.name === "Earth");
    const stats = document.getElementById('simulationStats');
    if (earth) {
        const velocity = earth.velocity.magnitude();
        stats.innerHTML = `
            Earth Velocity: ${(velocity).toFixed(2)} AU/year<br>
            Distance from Sun: ${earth.position.magnitude().toFixed(3)} AU<br>
            Relativistic Factor: ${(1 / Math.sqrt(1 - (velocity * velocity) / (c * c))).toFixed(6)}
        `;
    }
}

// Event listeners
document.getElementById('toggleSimulation').addEventListener('click', () => {
    isRunning = !isRunning;
    if (isRunning) {
        animate();
    } else {
        cancelAnimationFrame(animationFrameId);
    }
});

document.getElementById('resetSimulation').addEventListener('click', () => {
    // Reset all bodies to initial positions
    bodies.length = 0;
    bodies.push(
        new CelestialBody("Sun", 1.0, 0.04, new Vector3D(0, 0, 0), new Vector3D(0, 0, 0), "#FFD700", true),
        new CelestialBody("Mercury", 1.651e-7, 0.015, new Vector3D(0.387, 0, 0), new Vector3D(0, 10.06, 0), "#A0522D"),
        new CelestialBody("Venus", 2.447e-6, 0.02, new Vector3D(0.723, 0, 0), new Vector3D(0, 7.39, 0), "#DEB887"),
        new CelestialBody("Earth", 3.003e-6, 0.02, new Vector3D(1, 0, 0), new Vector3D(0, 6.283, 0), "#4169E1"),
        new CelestialBody("Mars", 3.213e-7, 0.018, new Vector3D(1.524, 0, 0), new Vector3D(0, 5.09, 0), "#CD5C5C")
    );
    addMoonsAndAsteroids();
    
    // Reset light rays
    lightRays.length = 0;
    initLightRays();
});

// Handle window resize
window.addEventListener('resize', () => {
    canvas.width = window.innerWidth * 0.8;
    canvas.height = window.innerHeight * 0.8;
}); 
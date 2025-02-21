package qwe;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;



public class SolarSystemTopDown extends JFrame {
    public SolarSystemTopDown() {
        setTitle("Top-Down 2D Solar System Simulation with 3D Orbits, General Relativity, and Interactive Space-Time Bending");
        setSize(800, 600);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setLocationRelativeTo(null);
        SimulationPanel2D panel = new SimulationPanel2D();
        add(panel);
        setVisible(true);
        panel.start();
    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> new SolarSystemTopDown());
    }
}

class SimulationPanel2D extends JPanel implements KeyListener, Runnable {
    private final double G = 39.47841760435743;
    private final double c = 63239.7263;
    private final double dt = 0.001;

    private java.util.List<CelestialBody> bodies;
    private java.util.List<LightRay> lightRays;
    private Thread simulationThread;
    private boolean running = false;
    private double scale = 100;  
    private int panelWidth, panelHeight;

    public SimulationPanel2D() {
        setBackground(Color.BLACK);
        setFocusable(true);
        requestFocusInWindow();
        addKeyListener(this);
        initBodies();
        initLightRays();
    }

    private void initBodies() {
        bodies = new ArrayList<>();
        bodies.add(new CelestialBody("Sun", 1.0,
                new Vector3D(0, 0, 0), new Vector3D(0.01, 0, 0), 8, true, Color.YELLOW));
        
        bodies.add(new CelestialBody("Mercury", 1.651e-7,
                new Vector3D(0.39, 0.05, 0), new Vector3D(0, 10.06, 2), 3, false, Color.LIGHT_GRAY));
        
        bodies.add(new CelestialBody("Venus", 2.447e-6,
                new Vector3D(0.723, -0.03, 0.1), new Vector3D(0, 7.39, -1), 4, false, Color.PINK));
        
        bodies.add(new CelestialBody("Earth", 3.003e-6,
                new Vector3D(1, 0, 0), new Vector3D(0, 6.283, 0), 4, false, Color.BLUE));
        
        bodies.add(new CelestialBody("Mars", 3.213e-7,
                new Vector3D(1.524, 0.1, -0.2), new Vector3D(0, 5.09, 0.5), 3, false, Color.RED));
        
        bodies.add(new CelestialBody("Jupiter", 9.545e-4,
                new Vector3D(5.203, -0.2, 0.3), new Vector3D(0, 2.76, -0.3), 6, false, Color.ORANGE));
        
        bodies.add(new CelestialBody("Saturn", 2.858e-4,
                new Vector3D(9.537, 0.15, 0.4), new Vector3D(0, 2.04, 0.2), 5, false, Color.LIGHT_GRAY));


        CelestialBody earth = getBodyByName("Earth");
        if (earth != null) {
            double distance = 0.05;
            double orbitalSpeed = Math.sqrt(G * earth.mass / distance);
            CelestialBody earthMoon = new CelestialBody("EarthMoon", 3.7e-8,
                    earth.position.add(new Vector3D(distance, 0, 0)),
                    earth.velocity.add(new Vector3D(0, orbitalSpeed, 0)),
                    2, false, Color.LIGHT_GRAY);
            bodies.add(earthMoon);
        }
        CelestialBody mars = getBodyByName("Mars");
        if (mars != null) {
            double distance = 0.03;
            double orbitalSpeed = Math.sqrt(G * mars.mass / distance);
            CelestialBody marsMoon1 = new CelestialBody("MarsMoon1", 1e-8,
                    mars.position.add(new Vector3D(distance, 0, 0)),
                    mars.velocity.add(new Vector3D(0, orbitalSpeed, 0)),
                    1, false, Color.GRAY);
            CelestialBody marsMoon2 = new CelestialBody("MarsMoon2", 1e-8,
                    mars.position.add(new Vector3D(-distance, 0, 0)),
                    mars.velocity.add(new Vector3D(0, -orbitalSpeed, 0)),
                    1, false, Color.GRAY);
            bodies.add(marsMoon1);
            bodies.add(marsMoon2);
        }
        CelestialBody jupiter = getBodyByName("Jupiter");
        if (jupiter != null) {
            double distance = 0.1;
            double orbitalSpeed = Math.sqrt(G * jupiter.mass / distance);
            CelestialBody jupiterMoon1 = new CelestialBody("JupiterMoon1", 5e-8,
                    jupiter.position.add(new Vector3D(distance, 0, 0)),
                    jupiter.velocity.add(new Vector3D(0, orbitalSpeed, 0)),
                    2, false, Color.WHITE);
            CelestialBody jupiterMoon2 = new CelestialBody("JupiterMoon2", 5e-8,
                    jupiter.position.add(new Vector3D(0, 0, distance)),
                    jupiter.velocity.add(new Vector3D(-orbitalSpeed, 0, 0)),
                    2, false, Color.WHITE);
            CelestialBody jupiterMoon3 = new CelestialBody("JupiterMoon3", 5e-8,
                    jupiter.position.add(new Vector3D(-distance, 0, 0)),
                    jupiter.velocity.add(new Vector3D(0, -orbitalSpeed, 0)),
                    2, false, Color.WHITE);
            bodies.add(jupiterMoon1);
            bodies.add(jupiterMoon2);
            bodies.add(jupiterMoon3);
        }
        CelestialBody saturn = getBodyByName("Saturn");
        if (saturn != null) {
            double distance = 0.1;
            double orbitalSpeed = Math.sqrt(G * saturn.mass / distance);
            CelestialBody saturnMoon1 = new CelestialBody("SaturnMoon1", 5e-8,
                    saturn.position.add(new Vector3D(distance, 0, 0)),
                    saturn.velocity.add(new Vector3D(0, orbitalSpeed, 0)),
                    2, false, Color.LIGHT_GRAY);
            CelestialBody saturnMoon2 = new CelestialBody("SaturnMoon2", 5e-8,
                    saturn.position.add(new Vector3D(0, 0, distance)),
                    saturn.velocity.add(new Vector3D(-orbitalSpeed, 0, 0)),
                    2, false, Color.LIGHT_GRAY);
            CelestialBody saturnMoon3 = new CelestialBody("SaturnMoon3", 5e-8,
                    saturn.position.add(new Vector3D(-distance, 0, 0)),
                    saturn.velocity.add(new Vector3D(0, -orbitalSpeed, 0)),
                    2, false, Color.LIGHT_GRAY);
            bodies.add(saturnMoon1);
            bodies.add(saturnMoon2);
            bodies.add(saturnMoon3);
        }
        int numAsteroids = 100;
        CelestialBody sun = getBodyByName("Sun");
        if (sun != null) {
            Random rand = new Random();
            for (int i = 0; i < numAsteroids; i++) {
                double angle = rand.nextDouble() * 2 * Math.PI;
                double r = 2.2 + rand.nextDouble() * (3.2 - 2.2);
                double x = r * Math.cos(angle);
                double z = r * Math.sin(angle);
                double orbitalSpeed = Math.sqrt(G * sun.mass / r);
                double vx = -orbitalSpeed * Math.sin(angle);
                double vz = orbitalSpeed * Math.cos(angle);
                CelestialBody asteroid = new CelestialBody("Asteroid" + i, 1e-9,
                        new Vector3D(x, 0, z),
                        new Vector3D(vx, orbitalSpeed * 0.1, vz),
                        1, false, Color.DARK_GRAY);
                bodies.add(asteroid);
            }
        }
    }

    private CelestialBody getBodyByName(String name) {
        for (CelestialBody body : bodies) {
            if (body.name.equals(name)) {
                return body;
            }
        }
        return null;
    }

    private void initLightRays() {
        lightRays = new ArrayList<>();
        int numRays = 20;
        double startX = -5;
        double minZ = -3, maxZ = 3;
        for (int i = 0; i < numRays; i++) {
            double z = minZ + i * (maxZ - minZ) / (numRays - 1);
            LightRay ray = new LightRay(new Vector3D(startX, 0, z), new Vector3D(c, 0, 0));
            lightRays.add(ray);
        }
    }

    public void start() {
        running = true;
        simulationThread = new Thread(this);
        simulationThread.start();
    }

    @Override
    public void run() {
        while (running) {
            updateSimulationRK4();
            updateLightRaysRK4();
            repaint();
            try {
                Thread.sleep(16);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
    }

    private void updateSimulationRK4() {
        int n = bodies.size();
        Vector3D[] pos = new Vector3D[n];
        Vector3D[] vel = new Vector3D[n];
        for (int i = 0; i < n; i++) {
            pos[i] = bodies.get(i).position;
            vel[i] = bodies.get(i).velocity;
        }
        
        Vector3D[] dp1 = new Vector3D[n];
        Vector3D[] dv1 = new Vector3D[n];
        for (int i = 0; i < n; i++) {
            Vector3D a = computeAcceleration(i, pos, vel);
            dp1[i] = vel[i].multiply(dt);
            dv1[i] = a.multiply(dt);
        }
        
        Vector3D[] pos2 = new Vector3D[n];
        Vector3D[] vel2 = new Vector3D[n];
        for (int i = 0; i < n; i++) {
            pos2[i] = pos[i].add(dp1[i].multiply(0.5));
            vel2[i] = vel[i].add(dv1[i].multiply(0.5));
        }
        Vector3D[] dp2 = new Vector3D[n];
        Vector3D[] dv2 = new Vector3D[n];
        for (int i = 0; i < n; i++) {
            Vector3D a = computeAcceleration(i, pos2, vel2);
            dp2[i] = vel2[i].multiply(dt);
            dv2[i] = a.multiply(dt);
        }
        
        Vector3D[] pos3 = new Vector3D[n];
        Vector3D[] vel3 = new Vector3D[n];
        for (int i = 0; i < n; i++) {
            pos3[i] = pos[i].add(dp2[i].multiply(0.5));
            vel3[i] = vel[i].add(dv2[i].multiply(0.5));
        }
        Vector3D[] dp3 = new Vector3D[n];
        Vector3D[] dv3 = new Vector3D[n];
        for (int i = 0; i < n; i++) {
            Vector3D a = computeAcceleration(i, pos3, vel3);
            dp3[i] = vel3[i].multiply(dt);
            dv3[i] = a.multiply(dt);
        }
        
        Vector3D[] pos4 = new Vector3D[n];
        Vector3D[] vel4 = new Vector3D[n];
        for (int i = 0; i < n; i++) {
            pos4[i] = pos[i].add(dp3[i]);
            vel4[i] = vel[i].add(dv3[i]);
        }
        Vector3D[] dp4 = new Vector3D[n];
        Vector3D[] dv4 = new Vector3D[n];
        for (int i = 0; i < n; i++) {
            Vector3D a = computeAcceleration(i, pos4, vel4);
            dp4[i] = vel4[i].multiply(dt);
            dv4[i] = a.multiply(dt);
        }
        
        for (int i = 0; i < n; i++) {
            Vector3D dp = dp1[i].add(dp2[i].multiply(2)).add(dp3[i].multiply(2)).add(dp4[i]).multiply(1.0 / 6.0);
            Vector3D dv = dv1[i].add(dv2[i].multiply(2)).add(dv3[i].multiply(2)).add(dv4[i]).multiply(1.0 / 6.0);
            bodies.get(i).position = bodies.get(i).position.add(dp);
            bodies.get(i).velocity = bodies.get(i).velocity.add(dv);
        }
    }

    private Vector3D computeAcceleration(int i, Vector3D[] pos, Vector3D[] vel) {
        Vector3D acc = new Vector3D(0, 0, 0);
        CelestialBody bi = bodies.get(i);
        for (int j = 0; j < pos.length; j++) {
            if (i == j)
                continue;
            Vector3D rVec = pos[i].subtract(pos[j]);
            double r = rVec.magnitude();
            if (r < 1e-6)
                continue;
            double factor = G * bodies.get(j).mass / (r * r * r);
            if (bodies.get(j).isSun && !bi.isSun) {
                Vector3D relVel = vel[i].subtract(vel[j]);
                double l = rVec.cross(relVel).magnitude();
                double correction = 1 + 3 * (l * l) / (r * r * c * c);
                factor *= correction;
            }
            acc = acc.subtract(rVec.multiply(factor));
        }
        return acc;
    }

    private void updateLightRaysRK4() {
        for (LightRay ray : lightRays) {
            rk4LightRay(ray);
            if (ray.position.x > 8) {
                ray.reset(-5, ray.position.z, c);
            }
        }
    }

    private void rk4LightRay(LightRay ray) {
        Vector3D p0 = ray.position;
        Vector3D v0 = ray.velocity;
        Vector3D a1 = computeLightAcceleration(ray);
        Vector3D dp1 = v0.multiply(dt);
        Vector3D dv1 = a1.multiply(dt);
        
        Vector3D p2 = p0.add(dp1.multiply(0.5));
        Vector3D v2 = v0.add(dv1.multiply(0.5));
        LightRay tempRay = new LightRay(p2, v2);
        Vector3D a2 = computeLightAcceleration(tempRay);
        Vector3D dp2 = v2.multiply(dt);
        Vector3D dv2 = a2.multiply(dt);
        
        Vector3D p3 = p0.add(dp2.multiply(0.5));
        Vector3D v3 = v0.add(dv2.multiply(0.5));
        tempRay = new LightRay(p3, v3);
        Vector3D a3 = computeLightAcceleration(tempRay);
        Vector3D dp3 = v3.multiply(dt);
        Vector3D dv3 = a3.multiply(dt);
        
        Vector3D p4 = p0.add(dp3);
        Vector3D v4 = v0.add(dv3);
        tempRay = new LightRay(p4, v4);
        Vector3D a4 = computeLightAcceleration(tempRay);
        Vector3D dp4 = v4.multiply(dt);
        Vector3D dv4 = a4.multiply(dt);
        
        Vector3D dp = dp1.add(dp2.multiply(2)).add(dp3.multiply(2)).add(dp4).multiply(1.0 / 6.0);
        Vector3D dv = dv1.add(dv2.multiply(2)).add(dv3.multiply(2)).add(dv4).multiply(1.0 / 6.0);
        
        ray.position = ray.position.add(dp);
        ray.velocity = ray.velocity.add(dv).normalize().multiply(c);
        ray.path.add(new Vector3D(ray.position.x, ray.position.y, ray.position.z));
        if (ray.path.size() > 500) {
            ray.path.remove(0);
        }
    }

    private Vector3D computeLightAcceleration(LightRay ray) {
        Vector3D a = new Vector3D(0, 0, 0);
        Vector3D n = ray.velocity.normalize();
        for (CelestialBody body : bodies) {
            Vector3D d = ray.position.subtract(body.position);
            double r = d.magnitude();
            if (r < 1e-6)
                continue;
            Vector3D dParallel = n.multiply(d.dot(n));
            Vector3D dPerp = d.subtract(dParallel);
            if (dPerp.magnitude() < 1e-6)
                continue;
            Vector3D nPerp = dPerp.normalize();
            Vector3D ai = nPerp.multiply((4 * G * body.mass) / (r * r * c));
            a = a.add(ai);
        }
        return a;
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        panelWidth = getWidth();
        panelHeight = getHeight();
        Graphics2D g2d = (Graphics2D) g;
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        drawGrid(g2d, bodies);

        for (CelestialBody body : bodies) {
            double screenX = body.position.x * scale + panelWidth / 2;
            double screenY = body.position.z * scale + panelHeight / 2;
            double drawRadius = body.radius;
            g2d.setColor(body.color);
            g2d.fillOval((int) (screenX - drawRadius), (int) (screenY - drawRadius),
                    (int) (2 * drawRadius), (int) (2 * drawRadius));
        }

        g2d.setColor(Color.CYAN);
        for (LightRay ray : lightRays) {
            if (ray.path.size() > 1) {
                for (int i = 0; i < ray.path.size() - 1; i++) {
                    Vector3D p1 = ray.path.get(i);
                    Vector3D p2 = ray.path.get(i + 1);
                    double screenX1 = p1.x * scale + panelWidth / 2;
                    double screenY1 = p1.z * scale + panelHeight / 2;
                    double screenX2 = p2.x * scale + panelWidth / 2;
                    double screenY2 = p2.z * scale + panelHeight / 2;
                    g2d.drawLine((int) screenX1, (int) screenY1, (int) screenX2, (int) screenY2);
                }
            }
        }

        g2d.setColor(Color.WHITE);
        g2d.drawString("Zoom: Up/Down Arrow Keys", 10, 20);
    }

    private void drawGrid(Graphics2D g2d, java.util.List<CelestialBody> bodies) {
        g2d.setColor(new Color(50, 50, 50));
        double gridSpacing = 0.2;
        double simLeft = -panelWidth / (2.0 * scale);
        double simRight = panelWidth / (2.0 * scale);
        double simTop = -panelHeight / (2.0 * scale);
        double simBottom = panelHeight / (2.0 * scale);

        for (double x = Math.floor(simLeft / gridSpacing) * gridSpacing; x <= simRight; x += gridSpacing) {
            java.util.List<Point> polyline = new ArrayList<>();
            for (double z = simTop; z <= simBottom; z += gridSpacing / 10) {
                Vector3D original = new Vector3D(x, 0, z);
                Vector3D bent = bendPoint(original, bodies);
                int screenX = (int) (bent.x * scale + panelWidth / 2);
                int screenY = (int) (bent.z * scale + panelHeight / 2);
                polyline.add(new Point(screenX, screenY));
            }
            for (int i = 0; i < polyline.size() - 1; i++) {
                Point p1 = polyline.get(i);
                Point p2 = polyline.get(i + 1);
                g2d.drawLine(p1.x, p1.y, p2.x, p2.y);
            }
        }

        for (double z = Math.floor(simTop / gridSpacing) * gridSpacing; z <= simBottom; z += gridSpacing) {
            java.util.List<Point> polyline = new ArrayList<>();
            for (double x = simLeft; x <= simRight; x += gridSpacing / 10) {
                Vector3D original = new Vector3D(x, 0, z);
                Vector3D bent = bendPoint(original, bodies);
                int screenX = (int) (bent.x * scale + panelWidth / 2);
                int screenY = (int) (bent.z * scale + panelHeight / 2);
                polyline.add(new Point(screenX, screenY));
            }
            for (int i = 0; i < polyline.size() - 1; i++) {
                Point p1 = polyline.get(i);
                Point p2 = polyline.get(i + 1);
                g2d.drawLine(p1.x, p1.y, p2.x, p2.y);
            }
        }
    }

    private Vector3D bendPoint(Vector3D p, java.util.List<CelestialBody> bodies) {
        Vector3D offset = new Vector3D(0, 0, 0);
        double k = 0.1;
        for (CelestialBody body : bodies) {
            Vector3D d = p.subtract(body.position);
            double r = d.magnitude();
            double minR = 0.2;
            if (r < minR)
                r = minR;
            offset = offset.add(d.normalize().multiply(k * body.mass / r));
        }
        return p.subtract(offset);
    }

    @Override
    public void keyPressed(KeyEvent e) {
        int key = e.getKeyCode();
        if (key == KeyEvent.VK_UP) {
            scale += 10;
        } else if (key == KeyEvent.VK_DOWN) {
            scale = Math.max(10, scale - 10);
        }
    }

    @Override
    public void keyReleased(KeyEvent e) { }
    @Override
    public void keyTyped(KeyEvent e) { }
}

class CelestialBody {
    String name;
    double mass;
    Vector3D position;
    Vector3D velocity;
    double radius;
    boolean isSun;
    Color color;

    public CelestialBody(String name, double mass, Vector3D position, Vector3D velocity,
                           double radius, boolean isSun, Color color) {
        this.name = name;
        this.mass = mass;
        this.position = position;
        this.velocity = velocity;
        this.radius = radius;
        this.isSun = isSun;
        this.color = color;
    }
}

class Vector3D {
    double x, y, z;
    public Vector3D(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
    public Vector3D add(Vector3D v) {
        return new Vector3D(x + v.x, y + v.y, z + v.z);
    }
    public Vector3D subtract(Vector3D v) {
        return new Vector3D(x - v.x, y - v.y, z - v.z);
    }
    public Vector3D multiply(double scalar) {
        return new Vector3D(x * scalar, y * scalar, z * scalar);
    }
    public double dot(Vector3D v) {
        return x * v.x + y * v.y + z * v.z;
    }
    public Vector3D cross(Vector3D v) {
        return new Vector3D(
            y * v.z - z * v.y,
            z * v.x - x * v.z,
            x * v.y - y * v.x
        );
    }
    public double magnitude() {
        return Math.sqrt(x * x + y * y + z * z);
    }
    public Vector3D normalize() {
        double m = magnitude();
        if (m == 0)
            return new Vector3D(0, 0, 0);
        return new Vector3D(x / m, y / m, z / m);
    }
}

class LightRay {
    Vector3D position;
    Vector3D velocity;
    ArrayList<Vector3D> path;

    public LightRay(Vector3D position, Vector3D velocity) {
        this.position = position;
        this.velocity = velocity.normalize().multiply(velocity.magnitude());
        this.path = new ArrayList<>();
        this.path.add(new Vector3D(position.x, position.y, position.z));
    }

    public void reset(double startX, double z, double c) {
        this.position = new Vector3D(startX, 0, z);
        this.velocity = new Vector3D(c, 0, 0);
        this.path.clear();
        this.path.add(new Vector3D(this.position.x, this.position.y, this.position.z));
    }
}
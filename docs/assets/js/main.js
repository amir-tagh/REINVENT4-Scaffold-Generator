"use strict";

/* =========================================
   THEME MANAGEMENT
   ========================================= */

const root = document.documentElement;
const themeToggle = document.getElementById("theme-toggle");

function setTheme(theme) {
    root.setAttribute("data-theme", theme);

    try {
        localStorage.setItem("reinvent-theme", theme);
    } catch (error) {
        // Continue when browser storage is unavailable.
    }

    if (themeToggle) {
        themeToggle.textContent =
            theme === "dark" ? "\u2600" : "\u263E";

        themeToggle.setAttribute(
            "aria-label",
            theme === "dark"
                ? "Switch to light mode"
                : "Switch to dark mode"
        );
    }
}

let savedTheme = "dark";

try {
    savedTheme = localStorage.getItem("reinvent-theme") || "dark";
} catch (error) {
    // Use the default theme.
}

setTheme(savedTheme === "light" ? "light" : "dark");

if (themeToggle) {
    themeToggle.addEventListener("click", () => {
        const current = root.getAttribute("data-theme");
        setTheme(current === "dark" ? "light" : "dark");
    });
}


/* =========================================
   MOBILE NAVIGATION
   ========================================= */

const menuToggle = document.getElementById("menu-toggle");
const navLinks = document.getElementById("nav-links");

if (menuToggle && navLinks) {
    menuToggle.addEventListener("click", () => {
        const isOpen = navLinks.classList.toggle("open");

        menuToggle.setAttribute(
            "aria-expanded",
            String(isOpen)
        );
    });

    navLinks.querySelectorAll("a").forEach(link => {
        link.addEventListener("click", () => {
            navLinks.classList.remove("open");
            menuToggle.setAttribute("aria-expanded", "false");
        });
    });

    document.addEventListener("keydown", event => {
        if (event.key === "Escape") {
            navLinks.classList.remove("open");
            menuToggle.setAttribute("aria-expanded", "false");
        }
    });
}


/* =========================================
   CMOS V2 OPTIMIZATION COMPONENTS
   ========================================= */

const components = [
    { name: "Drug-likeness (QED)", color: "#38bdf8" },
    { name: "Synthetic Accessibility", color: "#818cf8" },
    { name: "Synthetic Complexity", color: "#a78bfa" },
    { name: "Parent Similarity", color: "#e879f9" },
    { name: "Structural Diversity", color: "#2dd4bf" },
    { name: "Solubility", color: "#22c55e" },
    { name: "CNS Desirability", color: "#facc15" },
    { name: "BBB Desirability", color: "#fb923c" },
    { name: "ADME Desirability", color: "#f472b6" },
    { name: "Toxicity Desirability", color: "#f87171" }
];


/* =========================================
   CMOS V2 OPTIMIZATION PROFILES
   ========================================= */


const componentKeys = [
    "QED",
    "SA",
    "SCScore",
    "ParentSimilarity",
    "Diversity",
    "Solubility",
    "CNS",
    "BBB",
    "ADME",
    "Toxicity"
];

const profileMetadata = {
    cns: {
        label: "CNS Optimization",
        description: "Default CNS-oriented molecular optimization."
    },
    general: {
        label: "General Optimization",
        description: "General optimization without CNS or BBB objectives."
    },
    oral: {
        label: "Oral Drug Discovery",
        description: "Emphasizes drug-likeness, solubility and ADME."
    },
    experimental: {
        label: "Experimental Prioritization",
        description: "Emphasizes diversity, solubility and synthetic accessibility."
    },
    cns_original: {
        label: "CNS Original",
        description: "Original CNS weighting configuration."
    },
    cns_toxicity: {
        label: "CNS Toxicity",
        description: "Increases toxicity desirability weighting relative to CNS Original."
    },
    cns_safety: {
        label: "CNS Safety",
        description: "Safety-oriented configuration identical to the default CNS profile."
    },
    cns_diversity: {
        label: "CNS Diversity",
        description: "Emphasizes structural diversity during molecular selection."
    },
    cns_balanced: {
        label: "CNS Balanced",
        description: "Balances ADME and toxicity desirability objectives."
    }
};

const profiles = {};

const profileSelect = document.getElementById("profile-select");
const profileDescription = document.getElementById("profile-description");
const weightDonut = document.getElementById("weight-donut");
const weightLegend = document.getElementById("weight-legend");
const weightBars = document.getElementById("weight-bars");

function renderProfile(profileName) {
    const profile = profiles[profileName];

    if (!profile || !weightDonut || !weightLegend || !weightBars) {
        return;
    }

    const weights = profile.weights;
    const total = weights.reduce((sum, weight) => sum + weight, 0);

    if (profileDescription) {
        profileDescription.textContent = profile.description;
    }

    /* Donut chart */

    let start = 0;
    const segments = [];

    weights.forEach((weight, index) => {
        if (weight === 0) return;

        const end = start + (weight / total) * 100;

        segments.push(
            `${components[index].color} ${start}% ${end}%`
        );

        start = end;
    });

    weightDonut.style.background =
        `conic-gradient(${segments.join(",")})`;

    weightDonut.querySelector("strong").textContent =
        `${total}%`;

    weightDonut.setAttribute(
        "aria-label",
        `Optimization weights for ${profileName}`
    );

    /* Chart legend */

    weightLegend.innerHTML = "";

    components.forEach((component, index) => {
        const item = document.createElement("div");
        item.className = "legend-item";

        const color = document.createElement("span");
        color.className = "legend-color";
        color.style.background = component.color;

        const label = document.createElement("span");
        label.textContent =
            `${component.name}: ${weights[index]}%`;

        item.append(color, label);
        weightLegend.appendChild(item);
    });

    /* Horizontal weight bars */

    weightBars.innerHTML = "";

    components.forEach((component, index) => {
        const weight = weights[index];

        const row = document.createElement("div");
        row.className = "weight-row";

        const label = document.createElement("div");
        label.className = "weight-label";

        const name = document.createElement("span");
        name.textContent = component.name;

        const value = document.createElement("span");
        value.className = "weight-value";
        value.textContent = `${weight}%`;

        label.append(name, value);

        const track = document.createElement("div");
        track.className = "weight-track";

        const fill = document.createElement("div");
        fill.className = "weight-fill";
        fill.style.background = component.color;
        fill.style.width = `${weight}%`;

        track.setAttribute("role", "meter");
        track.setAttribute("aria-label", component.name);
        track.setAttribute("aria-valuemin", "0");
        track.setAttribute("aria-valuemax", "100");
        track.setAttribute("aria-valuenow", String(weight));

        track.appendChild(fill);
        row.append(label, track);
        weightBars.appendChild(row);
    });
}

async function initializeProfiles() {
    if (!profileSelect) return;

    profileSelect.disabled = true;

    try {
        const response = await fetch(
            "assets/data/optimization_profiles.json",
            { cache: "no-cache" }
        );

        if (!response.ok) {
            throw new Error("Configuration could not be loaded.");
        }

        const configuration = await response.json();

        profileSelect.innerHTML = "";

        Object.entries(configuration).forEach(([name, weights]) => {
            const metadata = profileMetadata[name] || {
                label: name,
                description: "Custom optimization profile."
            };

            profiles[name] = {
                description: metadata.description,
                weights: componentKeys.map(
                    key => Number((weights[key] * 100).toFixed(6))
                )
            };

            const option = document.createElement("option");
            option.value = name;
            option.textContent = metadata.label;

            profileSelect.appendChild(option);
        });

        if (profiles.cns) {
            profileSelect.value = "cns";
        }

        profileSelect.disabled = false;

        profileSelect.addEventListener("change", event => {
            renderProfile(event.target.value);
        });

        renderProfile(profileSelect.value);

        console.info(
            "Loaded " + Object.keys(profiles).length +
            " optimization profiles."
        );

    } catch (error) {
        console.error(error);

        if (profileDescription) {
            profileDescription.textContent =
                "Unable to load optimization profiles.";
        }
    }
}

initializeProfiles();


/* =========================================
   COPY INSTALLATION COMMANDS
   ========================================= */

document.querySelectorAll(".copy-button").forEach(button => {
    button.addEventListener("click", async () => {
        const target = document.getElementById(
            button.dataset.copy
        );

        if (!target) return;

        try {
            await navigator.clipboard.writeText(
                target.textContent
            );

            button.textContent = "Copied!";

        } catch (error) {
            button.textContent = "Copy manually";
        }

        setTimeout(() => {
            button.textContent = "Copy";
        }, 1800);
    });
});


/* =========================================
   ANIMATED MOLECULAR NETWORK
   Decorative scientific visualization
   ========================================= */

const canvas = document.getElementById("molecular-canvas");

if (canvas && canvas.getContext) {
    const ctx = canvas.getContext("2d");

    const reducedMotion = window.matchMedia(
        "(prefers-reduced-motion: reduce)"
    );

    let width = 0;
    let height = 0;
    let scale = 1;
    let frameId = null;

    const rings = [
        { x: -1.5, y: 0.0, type: 0 },
        { x: -0.5, y: -0.87, type: 0 },
        { x: 0.5, y: 0.0, type: 1 },
        { x: 1.5, y: -0.87, type: 0 },
        { x: 1.5, y: 0.87, type: 1 },
        { x: -0.5, y: 0.87, type: 0 }
    ];

    const connections = [
        [0, 1],
        [0, 5],
        [1, 2],
        [2, 3],
        [2, 4],
        [2, 5]
    ];

    function resizeCanvas() {
        const rect = canvas.getBoundingClientRect();
        const ratio = Math.min(window.devicePixelRatio || 1, 2);

        width = rect.width;
        height = rect.height;

        canvas.width = width * ratio;
        canvas.height = height * ratio;

        ctx.setTransform(ratio, 0, 0, ratio, 0, 0);

        scale = Math.min(width / 7, height / 5);
    }

    function hexagon(x, y, radius, color) {
        ctx.beginPath();

        for (let i = 0; i < 6; i++) {
            const angle = Math.PI / 3 * i - Math.PI / 6;
            const px = x + radius * Math.cos(angle);
            const py = y + radius * Math.sin(angle);

            if (i === 0) {
                ctx.moveTo(px, py);
            } else {
                ctx.lineTo(px, py);
            }
        }

        ctx.closePath();
        ctx.fillStyle = color + "12";
        ctx.fill();

        ctx.strokeStyle = color;
        ctx.lineWidth = 1.6;
        ctx.stroke();

        /* Inner bond lines */

        for (let i = 0; i < 6; i += 2) {
            const a1 = Math.PI / 3 * i - Math.PI / 6;
            const a2 = Math.PI / 3 * (i + 1) - Math.PI / 6;

            ctx.beginPath();

            ctx.moveTo(
                x + radius * 0.73 * Math.cos(a1),
                y + radius * 0.73 * Math.sin(a1)
            );

            ctx.lineTo(
                x + radius * 0.73 * Math.cos(a2),
                y + radius * 0.73 * Math.sin(a2)
            );

            ctx.lineWidth = 1.2;
            ctx.stroke();
        }
    }

    function draw(time) {
        ctx.clearRect(0, 0, width, height);

        const light = root.getAttribute("data-theme") === "light";
        const blue = light ? "#0369a1" : "#38bdf8";
        const teal = light ? "#0f766e" : "#2dd4bf";

        const t = reducedMotion.matches ? 0 : time * 0.001;
        const radius = scale * 0.67;

        const positions = rings.map((ring, index) => ({
            x: width / 2 + ring.x * scale * 1.65,
            y: height / 2 + ring.y * scale * 1.65
                + Math.sin(t + index * 0.7) * 5
        }));

        /* Background glow */

        const glow = ctx.createRadialGradient(
            width / 2, height / 2, 0,
            width / 2, height / 2,
            Math.max(width, height) * 0.5
        );

        glow.addColorStop(0, blue + "15");
        glow.addColorStop(1, blue + "00");

        ctx.fillStyle = glow;
        ctx.fillRect(0, 0, width, height);

        /* Network connections */

        connections.forEach(([a, b]) => {
            const p1 = positions[a];
            const p2 = positions[b];

            ctx.beginPath();
            ctx.moveTo(p1.x, p1.y);
            ctx.lineTo(p2.x, p2.y);

            ctx.strokeStyle = blue + "55";
            ctx.lineWidth = 1.2;
            ctx.setLineDash([4, 5]);
            ctx.stroke();
        });

        ctx.setLineDash([]);

        /* Molecular rings */

        positions.forEach((position, index) => {
            const color = rings[index].type === 1
                ? teal
                : blue;

            hexagon(
                position.x,
                position.y,
                radius,
                color
            );

            const angle = Math.PI / 3 * index - Math.PI / 6;

            const atomX =
                position.x + radius * Math.cos(angle);

            const atomY =
                position.y + radius * Math.sin(angle);

            ctx.beginPath();
            ctx.arc(atomX, atomY, 3.5, 0, Math.PI * 2);
            ctx.fillStyle = color;
            ctx.fill();
        });

        /* Animated signal particles */

        connections.forEach(([a, b], index) => {
            const p1 = positions[a];
            const p2 = positions[b];

            const progress =
                (t * 0.22 + index * 0.17) % 1;

            const x = p1.x + (p2.x - p1.x) * progress;
            const y = p1.y + (p2.y - p1.y) * progress;

            ctx.beginPath();
            ctx.arc(x, y, 2.8, 0, Math.PI * 2);

            ctx.fillStyle = teal;
            ctx.shadowColor = teal;
            ctx.shadowBlur = 10;
            ctx.fill();
            ctx.shadowBlur = 0;
        });

        frameId = null;

        if (!reducedMotion.matches && !document.hidden) {
            frameId = requestAnimationFrame(draw);
        }
    }

    function restartAnimation() {
        if (frameId !== null) {
            cancelAnimationFrame(frameId);
        }

        frameId = null;
        resizeCanvas();
        draw(performance.now());
    }

    window.addEventListener("resize", restartAnimation);

    document.addEventListener("visibilitychange", () => {
        if (document.hidden) {
            if (frameId !== null) {
                cancelAnimationFrame(frameId);
                frameId = null;
            }
        } else {
            restartAnimation();
        }
    });

    reducedMotion.addEventListener("change", restartAnimation);

    if (themeToggle) {
        themeToggle.addEventListener("click", restartAnimation);
    }

    restartAnimation();
}


/* =========================================
   INITIALIZATION COMPLETE
   ========================================= */

console.info("Scientific platform website initialized.");

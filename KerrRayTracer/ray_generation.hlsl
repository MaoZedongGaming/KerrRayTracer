struct CameraConstants {
    float4 e0;
    float4 e1;
    float4 e2;
    float4 e3;
    
    float4 pos; 
    float4 vel;
    
    uint width;
    uint height;
    float fov;
    float a;
};

struct InitialPhoton {
    float r;
    float theta;
    float phi;
    float dr;
    
    float dtheta;
    float dphi;
    float xi;
    float eta;
    
    float sign_r;
    float sign_theta;
    uint state;
    float dlambda;
};

ConstantBuffer<CameraConstants> CameraData : register(b0);
RWStructuredBuffer<InitialPhoton> PhotonData : register(u0);

float g_tt(float r, float theta, float a) {
    float r2 = r * r;
    float cosTh = cos(theta);
    float sigma = r2 + a * a * cosTh * cosTh;
    return -(1.0f - (2.0f * r) / sigma);
}

float g_tphi(float r, float theta, float a) {
    float r2 = r * r;
    float sinTh = sin(theta);
    float cosTh = cos(theta);
    float sigma = r2 + a * a * cosTh * cosTh;
    return -(2.0f * a * r * sinTh * sinTh) / sigma;
}

float g_phiphi(float r, float theta, float a) {
    float r2 = r * r;
    float a2 = a * a;
    float sinTh = sin(theta);
    float sin2 = sinTh * sinTh;
    float cosTh = cos(theta);
    float sigma = r2 + a2 * cosTh * cosTh;
    return (r2 + a2 + (2.0f * a2 * r * sin2) / sigma) * sin2;
}

float g_thth(float r, float theta, float a) {
    float r2 = r * r;
    float cosTh = cos(theta);
    return r2 + a * a * cosTh * cosTh;
}

void evaluateInitialDerivatives(float r, float theta, float sign_r, float sign_theta, float xi, float eta, float a, out float dr, out float dtheta, out float dphi) {
    float r2 = r * r;
    float a2 = a * a;
    float cosTh = cos(theta);
    float sinTh = sin(theta);
    float sin2 = sinTh * sinTh;
    float delta = r2 - 2.0f * r + a2;
    float P = r2 + a2 - a * xi;
    float R = P * P - delta * (eta + (xi - a) * (xi - a));
    float Theta = eta + a2 * cosTh * cosTh - xi * xi / (sin2);
    
    dr = sign_r * sqrt(max(R, 0.0f));
    dtheta = sign_theta * sqrt(max(Theta, 0.0f));
    dphi = -(a - xi / sin2) + P * a / delta;
}


[numthreads(16, 16, 1)]
void CSMain(uint3 id : SV_DispatchThreadID) {
    if (id.x >= CameraData.width || id.y >= CameraData.height) {
        return;
    }
    uint photonIdx = id.x + id.y * CameraData.width;
    float screenX = (2.0f * ((float) id.x + 0.5f) / (float) CameraData.width - 1.0f) * tan(CameraData.fov / 2.0f) * ((float) CameraData.width / (float) CameraData.height);
    float screenY = (1.0f - 2.0f * ((float) id.y + 0.5f) / (float) CameraData.height) * tan(CameraData.fov / 2.0f);
    
    float a = CameraData.a;
    float r = CameraData.pos[1];
    float theta = CameraData.pos[2];
    float phi = CameraData.pos[3];
    
    float p_1 = 1.0f / sqrt(1.0f + screenX * screenX + screenY * screenY);
    float p_2 = screenY * p_1;
    float p_3 = screenX * p_1;
    
    float4 p = CameraData.e0 + -p_1 * CameraData.e1 + -p_2 * CameraData.e2 + p_3 * CameraData.e3;
    
    float g_t = g_tt(r, theta, a);
    float g_tp = g_tphi(r, theta, a);
    float g_pp = g_phiphi(r, theta, a);
    
    float cosTh = cos(theta);
    float sinTh = sin(theta);
    
    float E = -(p[0] * g_t + p[3] * g_tp);
    float L_z = p[3] * g_pp + p[0] * g_tp;
    float Q = p[2] * p[2] * g_t * g_t + cosTh * cosTh * (a * a * (-E * E) + (L_z * L_z) / (sinTh * sinTh));
    
    float xi = L_z / E;
    float eta = Q / (E * E);
    
    float sign_r = (p.y >= 0.0f) ? 1.0f : -1.0f;
    float sign_theta = (p.z >= 0.0f) ? 1.0f : -1.0f;

    float dr, dtheta, dphi;
    evaluateInitialDerivatives(r, theta, sign_r, sign_theta, xi, eta, a, dr, dtheta, dphi);

    InitialPhoton photon;
    photon.r = r;
    photon.theta = theta;
    photon.phi = phi;
    photon.dr = dr;
    photon.dtheta = dtheta;
    photon.dphi = dphi;
    photon.xi = xi;
    photon.eta = eta;
    photon.sign_r = sign_r;
    photon.sign_theta = sign_theta;
    photon.state = 0; // 0 = Active
    photon.dlambda = 0.1f; // Initial adaptive step guess

    PhotonData[photonIdx] = photon;
}
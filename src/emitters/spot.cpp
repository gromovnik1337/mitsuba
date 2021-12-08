/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/core/random.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/hw/basicshader.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{spot}{Spot light source}
 * \icon{emitter_spot}
 * \order{3}
 * \parameters{
 *     \parameter{toWorld}{\Transform\Or\Animation}{
 *        Specifies an optional sensor-to-world transformation.
 *        \default{none (i.e. sensor space $=$ world space)}
 *     }
 *     \parameter{intensity}{\Spectrum}{
 *         Specifies the maximum radiant intensity at the center
 *         in units of power per unit steradian.
 *         \default{1}
 *     }
 *     \parameter{cutoffAngle}{\Float}{Cutoff angle, beyond which the spot light is completely black \default{\code{20} degrees}}
 *     \parameter{beamWidth}{\Float}{Subtended angle of the central beam portion \default{\code{cutoffAngle}$\ \cdot\ \nicefrac 34$}}
 *     \parameter{texture}{\Texture}{
 *         An optional texture to be projected along the spot light
 *     }
 *     \parameter{samplingWeight}{\Float}{
 *         Specifies the relative amount of samples
 *         allocated to this emitter. \default{1}
 *     }
 * }
 *
 * This plugin provides a spot light with a linear falloff.
 * In its local coordinate system, the spot light is positioned at the origin
 * and points along the positive Z direction. It can be conveniently
 * reoriented using the \code{lookat} tag, e.g.:
 * \begin{xml}
 * <emitter type="spot">
 *     <transform name="toWorld">
 *         <!-- Orient the light so that points from (1, 1, 1) towards (1, 2, 1) -->
 *         <lookat origin="1, 1, 1" target="1, 2, 1"/>
 *     </transform>
 * </emitter>
 * \end{xml}
 *
 * The intensity linearly ramps up from \code{cutoffAngle}
 * to \code{beamWidth} (both specified in degrees), after which it remains at
 * the maximum value. A projection texture may optionally be supplied.
 */

ref<Random> g_random;

class SpotEmitter : public Emitter {
public:
    SpotEmitter(const Properties &props) : Emitter(props) {
        m_intensity = props.getSpectrum("intensity", Spectrum(1.0f));
        m_cutoffAngle = props.getFloat("cutoffAngle", 20);
        m_beamWidth = props.getFloat("beamWidth", m_cutoffAngle * 3.0f/4.0f);
        m_beamWidth = degToRad(m_beamWidth);
        m_cutoffAngle = degToRad(m_cutoffAngle);
        Assert(m_cutoffAngle >= m_beamWidth);
        m_type = EDeltaPosition;
        m_texture = new ConstantSpectrumTexture(
            props.getSpectrum("texture", Spectrum::getD65()));

        m_scaleX = props.getFloat("scaleX", 1.0f);
        m_scaleY = props.getFloat("scaleY", 1.0f);
        m_offX = props.getFloat("offX", 0.5f);
        m_offY = props.getFloat("offY", 0.5f);
        m_width = props.getInteger("width", 1);
        m_height = props.getInteger("height", 1);
        m_apertureRadius = props.getFloat("apertureRadius", 0.0f);
        m_focusDistance = props.getFloat("focusDistance", 1.0f);
        m_pixelGap = props.getFloat("pixelGap", 0.0f);

        /* Diffraction limit on resolution */
        m_diffLimit = props.getFloat("diffLimit", 0.0f);
        m_random = new Random();  // Requires extra (third) sample

        if (m_diffLimit < 0)
            Log(EError, "Diffraction limit on resolution (in degrees) "
                        "cannot be negative!");

        if (m_apertureRadius == 0) {
            Log(EWarn, "Can't have a zero aperture radius -- setting to %f", Epsilon / 10);
            m_apertureRadius = Epsilon / 10;
        }
        if (m_focusDistance == 0) {
            Log(EWarn, "Can't have a zero focus distance -- setting to %f", 1.0f);
            m_focusDistance = 1.0f;
        }
        if (m_pixelGap < 0) {
            Log(EWarn, "Pixel gap cannot be negative -- setting to %f", 0.0f);
            m_pixelGap = 0.0f;
        }
    }

    SpotEmitter(Stream *stream, InstanceManager *manager)
        : Emitter(stream, manager) {
        m_texture = static_cast<Texture *>(manager->getInstance(stream));
        m_intensity = Spectrum(stream);
        m_beamWidth = stream->readFloat();
        m_cutoffAngle = stream->readFloat();
        m_scaleX = stream->readFloat();
        m_scaleY = stream->readFloat();
        m_offX = stream->readFloat();
        m_offY = stream->readFloat();
        m_width = stream->readInt();
        m_height = stream->readInt();
        m_apertureRadius = stream->readFloat();
        m_focusDistance = stream->readFloat();
        m_diffLimit = stream->readFloat();
        m_pixelGap = stream->readFloat();
        configure();
    }

    void configure() {
        m_cosBeamWidth = std::cos(m_beamWidth);
        m_cosCutoffAngle = std::cos(m_cutoffAngle);
        m_uvFactor = std::tan(m_cutoffAngle);
        m_invTransitionWidth = 1.0f / (m_cutoffAngle - m_beamWidth);

        /* Use global copy of m_random to bypass const modifier of falloffCurve() function.
         * Such ad-hoc solution is potentially not thread safe
         * when using multiple projector light sources in one scene. */
        g_random = new Random(m_random);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Emitter::serialize(stream, manager);

        manager->serialize(stream, m_texture.get());
        m_intensity.serialize(stream);
        stream->writeFloat(m_beamWidth);
        stream->writeFloat(m_cutoffAngle);
        stream->writeFloat(m_scaleX);
        stream->writeFloat(m_scaleY);
        stream->writeFloat(m_offX);
        stream->writeFloat(m_offY);
        stream->writeInt(m_width);
        stream->writeInt(m_height);
        stream->writeFloat(m_apertureRadius);
        stream->writeFloat(m_focusDistance);
        stream->writeFloat(m_diffLimit);
        stream->writeFloat(m_pixelGap);
    }

    Spectrum samplePosition(PositionSamplingRecord &pRec, const Point2 &sample,
            const Point2 *extra) const {
        const Transform &trafo = m_worldTransform->eval(pRec.time);
        pRec.p = trafo.transformAffine(Point(0.0f));
        pRec.n = Normal(0.0f);
        pRec.pdf = 1.0f;
        pRec.measure = EDiscrete;
        return m_intensity * (4 * M_PI);
    }

    Spectrum evalPosition(const PositionSamplingRecord &pRec) const {
        return (pRec.measure == EDiscrete) ? (m_intensity * 4*M_PI) : Spectrum(0.0f);
    }

    Float pdfPosition(const PositionSamplingRecord &pRec) const {
        return (pRec.measure == EDiscrete) ? 1.0f : 0.0f;
    }

    Spectrum sampleDirection(DirectionSamplingRecord &dRec,
            PositionSamplingRecord &pRec,
            const Point2 &sample,
            const Point2 *extra) const {
        const Transform &trafo = m_worldTransform->eval(pRec.time);
        Vector d = warp::squareToUniformCone(m_cosCutoffAngle, sample);
        dRec.d = trafo(d);
        dRec.pdf = warp::squareToUniformConePdf(m_cosCutoffAngle);
        dRec.measure = ESolidAngle;
        return evalDirection(dRec, pRec)/dRec.pdf;
    }

    Float pdfDirection(const DirectionSamplingRecord &dRec,
            const PositionSamplingRecord &pRec) const {
        return (dRec.measure == ESolidAngle) ? warp::squareToUniformConePdf(m_cosCutoffAngle) : 0.0f;
    }

    Spectrum evalDirection(const DirectionSamplingRecord &dRec,
            const PositionSamplingRecord &pRec) const {
        const Transform &trafo = m_worldTransform->eval(pRec.time);
        return (dRec.measure == ESolidAngle) ?
            falloffCurve(trafo.inverse()(dRec.d)) * INV_FOURPI : Spectrum(0.0f);
    }

    Spectrum sampleRay(Ray &ray,
            const Point2 &spatialSample,
            const Point2 &directionalSample,
            Float time) const {
        const Transform &trafo = m_worldTransform->eval(time);

        Vector local = warp::squareToUniformCone(
            m_cosCutoffAngle, directionalSample);
        ray.setTime(time);
        ray.setOrigin(trafo.transformAffine(Point(0.0f)));
        ray.setDirection(trafo(local));
        Float dirPdf = warp::squareToUniformConePdf(m_cosCutoffAngle);
        return m_intensity * falloffCurve(local) / dirPdf;
    }

    Spectrum sampleDirect(DirectSamplingRecord &dRec, const Point2 &sample) const {
        const Transform &trafo = m_worldTransform->eval(dRec.time);

        Point2 tmp = warp::squareToUniformDiskConcentric(sample) * m_apertureRadius;
        Point apertureP(tmp.x, tmp.y, 0.0f);
        dRec.p = trafo.transformAffine(apertureP);
        dRec.pdf = 1.0f;
        dRec.measure = EDiscrete;
        dRec.uv = Point2(0.5f);
        dRec.d = dRec.p - dRec.ref;
        dRec.dist = dRec.d.length();
        Float invDist = 1.0f / dRec.dist;
        dRec.d *= invDist;
        dRec.n = Normal(0.0f);
        dRec.pdf = 1;
        dRec.measure = EDiscrete;

        Point refP = trafo.inverse().transformAffine(dRec.ref);
        Vector raydir = refP - apertureP;
        Vector focusP = Vector(apertureP) + raydir * (m_focusDistance / raydir.z);

        return m_intensity * falloffCurve(focusP) * (invDist * invDist);
    }

    inline Spectrum falloffCurve(const Vector &d) const {
        Spectrum result(1.0f);

        if (m_texture->getClass() != MTS_CLASS(ConstantSpectrumTexture)) {
            Intersection its;
            its.hasUVPartials = false;

            /* Project the point/direction onto image/mirror plane to get pixel coordinates for sampling */
            Point2 pix(m_scaleX * d.x / d.z + m_offX, m_scaleY * d.y / d.z + m_offY);

            /* Simulate diffraction limit on resolution */
            Point2 diffSample(g_random->nextFloat(), g_random->nextFloat());
            /* Do uniform sampling of Airy disk for simplicity. The diffraction limit effect is typically
             * on a scale of couple pixels so the difference is distributions (uniform vs actual diffraction profile)
             * will be indistiguishable (but the effect will become stronger so consider reducing the value of diffLimit) */
            Point2 diff = warp::squareToUniformDiskConcentric(diffSample);
            /* Convert diffraction limit from radians to pixels and add to the pixel sample */
            pix.x += diff.x * tan(degToRad(m_diffLimit)) * m_scaleX;
            pix.y += diff.y * tan(degToRad(m_diffLimit)) * m_scaleY;

            /* Evaluate texture only if the sample did not land in the gap between neighbour mirrors/pixels */
            Point2 frac(pix.x - floor(pix.x), pix.y - floor(pix.y));
            if (frac.x < m_pixelGap / 2 || frac.x > 1.0f - m_pixelGap / 2 ||
                frac.y < m_pixelGap / 2 || frac.y > 1.0f - m_pixelGap / 2)
                    return Spectrum(0.0f);

            its.uv = Point2(pix.x / m_width, pix.y / m_height);
            result = m_texture->eval(its);
        }
        return result;
    }

    Float pdfDirect(const DirectSamplingRecord &dRec) const {
        return dRec.measure == EDiscrete ? 1.0f : 0.0f;
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) && name == "texture") {
            m_texture = static_cast<Texture *>(child);
        } else {
            Emitter::addChild(name, child);
        }
    }

    AABB getAABB() const {
        //return m_worldTransform->getTranslationBounds();
        const Float r = m_apertureRadius;
        AABB bounds(Point(-r, -r, 0), Point(r, r, 0));
        return m_worldTransform->getSpatialBounds(bounds);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "SpotEmitter[" << std::endl
            << "  intensity = " << m_intensity.toString() << "," << std::endl
            << "  texture = " << m_texture.toString() << "," << std::endl
            << "  beamWidth = " << (m_beamWidth * 180/M_PI) << "," << std::endl
                << "  cutoffAngle = " << (m_cutoffAngle * 180/M_PI) << std::endl
                << "  scaleX = " << m_scaleX << "," << endl
                << "  scaleY = " << m_scaleY << "," << endl
                << "  offX = " << m_offX << "," << endl
                << "  offY = " << m_offY << "," << endl
                << "  width = " << m_width << "," << endl
                << "  height = " << m_height << "," << endl
                << "  apertureRadius = " << m_apertureRadius << "," << endl
                << "  focusDistance = " << m_focusDistance << "," << endl
                << "  diffLimit = " << m_diffLimit << "," << endl
                << "  pixelGap = " << m_pixelGap << "," << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
private:
    Spectrum m_intensity;
    ref<Texture> m_texture;
    ref<Random> m_random;
    Float m_beamWidth, m_cutoffAngle, m_uvFactor;
    Float m_cosBeamWidth, m_cosCutoffAngle, m_invTransitionWidth;
    Float m_scaleX, m_scaleY, m_offX, m_offY;
    Float m_width, m_height;
    Float m_apertureRadius;
    Float m_focusDistance;
    Float m_diffLimit;
    Float m_pixelGap;
};

// ================ Hardware shader implementation ================

class SpotEmitterShader : public Shader {
public:
    SpotEmitterShader(Renderer *renderer, Transform worldToEmitter,
        Float invTransitionWidth, Float cutoffAngle, Float cosCutoffAngle,
        Float cosBeamWidth, Float uvFactor, const Texture *texture)
        : Shader(renderer, EEmitterShader), m_worldToEmitter(worldToEmitter),
          m_invTransitionWidth(invTransitionWidth), m_cutoffAngle(cutoffAngle),
          m_cosCutoffAngle(cosCutoffAngle), m_cosBeamWidth(cosBeamWidth),
          m_uvFactor(uvFactor), m_texture(texture) {
        m_textureShader = renderer->registerShaderForResource(m_texture.get());
    }

    bool isComplete() const {
        return m_textureShader.get() != NULL;
    }

    void cleanup(Renderer *renderer) {
        renderer->unregisterShaderForResource(m_texture.get());
    }

    void putDependencies(std::vector<Shader *> &deps) {
        deps.push_back(m_textureShader.get());
    }

    void generateCode(std::ostringstream &oss, const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "uniform float " << evalName << "_invTransitionWidth;" << endl
            << "uniform float " << evalName << "_cutoffAngle;" << endl
            << "uniform float " << evalName << "_cosCutoffAngle;" << endl
            << "uniform float " << evalName << "_cosBeamWidth;" << endl
            << "uniform float " << evalName << "_uvFactor;" << endl
            << "uniform mat4 " << evalName << "_worldToEmitter;" << endl
            << "vec3 " << evalName << "_dir(vec3 wo) {" << endl
            << "    vec3 localDir = (" << evalName << "_worldToEmitter * vec4(wo, 0)).xyz;" << endl
            << "    float cosTheta = localDir.z;" << endl
            << "    if (cosTheta < " << evalName << "_cosCutoffAngle)" << endl
            << "        return vec3(0.0);" << endl
            << "    vec2 uv = 0.5 + 0.5 * (localDir.xy / (localDir.z * " << evalName << "_uvFactor));" << endl
            << "    vec3 color = " << depNames[0] << "(uv) * inv_fourpi;" << endl
            << "    if (cosTheta > " << evalName << "_cosBeamWidth)" << endl
            << "        return color;" << endl
            << "    return color * ((" << evalName << "_cutoffAngle - acos(cosTheta))" << endl
            << "           * " << evalName << "_invTransitionWidth);" << endl
            << "}" << endl;
    }

    void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
        parameterIDs.push_back(program->getParameterID(evalName + "_worldToEmitter", false));
        parameterIDs.push_back(program->getParameterID(evalName + "_invTransitionWidth", false));
        parameterIDs.push_back(program->getParameterID(evalName + "_cutoffAngle", false));
        parameterIDs.push_back(program->getParameterID(evalName + "_cosCutoffAngle", false));
        parameterIDs.push_back(program->getParameterID(evalName + "_cosBeamWidth", false));
        parameterIDs.push_back(program->getParameterID(evalName + "_uvFactor", false));
    }

    void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
        program->setParameter(parameterIDs[0], m_worldToEmitter);
        program->setParameter(parameterIDs[1], m_invTransitionWidth);
        program->setParameter(parameterIDs[2], m_cutoffAngle);
        program->setParameter(parameterIDs[3], m_cosCutoffAngle);
        program->setParameter(parameterIDs[4], m_cosBeamWidth);
        program->setParameter(parameterIDs[5], m_uvFactor);
    }

    MTS_DECLARE_CLASS()
private:
    Transform m_worldToEmitter;
    Float m_invTransitionWidth;
    Float m_cutoffAngle, m_cosCutoffAngle;
    Float m_cosBeamWidth, m_uvFactor;
    ref<const Texture> m_texture;
    ref<Shader> m_textureShader;
};

Shader *SpotEmitter::createShader(Renderer *renderer) const {
    const Transform &trafo = m_worldTransform->eval(0.0f);

    return new SpotEmitterShader(renderer, trafo.inverse(),
        m_invTransitionWidth, m_cutoffAngle, m_cosCutoffAngle,
        m_cosBeamWidth, m_uvFactor, m_texture.get());
}

MTS_IMPLEMENT_CLASS(SpotEmitterShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(SpotEmitter, false, Emitter)
MTS_EXPORT_PLUGIN(SpotEmitter, "Spot light");
MTS_NAMESPACE_END

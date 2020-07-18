#if defined(_WIN32)
#define _USE_MATH_DEFINES
#endif

#include <algorithm>

#include "avisynth.h"

enum KernelTypes
{
	AndrewsWave = 0,
	ElFallahFord,
	Gaussian,
	HubersMiniMax,
	Lorentzian,
	TukeyBiWeight,
	LinearDescent,
	Cosine,
	Flat,
	Inverse
};

enum ResTypes
{
	Mean = 0,
	Median,
	CWMedian,
	MultipleLinearRegression
};

class vsTBilateral : public GenericVideoFilter
{
	PClip ppclip_;
	int process[3];
	int diameter[3];
	double sDev[3], iDev[3], cs[3];
	bool d2_;
	int kerns_, kerni_, resType_;
	int y_, u_, v_;
	int pixel_max;
	double* spatialWeights[3];
	double* diffWeights[3];
	bool has_at_least_v8;

	template <typename PixelType>
	void ProcessFrameD2_Mean(const uint8_t* srcp_, const uint8_t* tp_, uint8_t* dstp_, int src_pitch, int tp_pitch, int dst_pitch, const int width, const int height, const int plane, IScriptEnvironment* env);
	template <typename PixelType>
	void ProcessFrameD2_MLR(const uint8_t* srcp_, const uint8_t* tp_, uint8_t* dstp_, int src_pitch, int tp_pitch, int dst_pitch, const int width, const int height, const int plane, IScriptEnvironment* env);
	template <typename PixelType>
	void ProcessFrameD2_Med(const uint8_t* srcp_, const uint8_t* tp_, uint8_t* dstp_, int src_pitch, int tp_pitch, int dst_pitch, const int width, const int height, const int plane, IScriptEnvironment* env);
	template <typename PixelType>
	void ProcessFrameD1_Mean(const uint8_t* srcp_, const uint8_t* tp_, uint8_t* dstp_, int src_pitch, int tp_pitch, int dst_pitch, const int width, const int height, const int plane, IScriptEnvironment* env);
	template <typename PixelType>
	void ProcessFrameD1_MLR(const uint8_t* srcp_, const uint8_t* tp_, uint8_t* dstp_, int src_pitch, int tp_pitch, int dst_pitch, const int width, const int height, const int plane, IScriptEnvironment* env);
	template <typename PixelType>
	void ProcessFrameD1_Med(const uint8_t* srcp_, const uint8_t* tp_, uint8_t* dstp_, int src_pitch, int tp_pitch, int dst_pitch, const int width, const int height, const int plane, IScriptEnvironment* env);

public:
	vsTBilateral(PClip _child, PClip ppclip, int diameterY, int diameterU, int diameterV, double sdevY, double sdevU, double sdevV, double idevY, double idevU,
		double idevV, double csY, double csU, double csV, bool d2, int kerns, int kerni, int restype, int y, int u, int v, IScriptEnvironment* env);
	PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);
	int __stdcall SetCacheHints(int cachehints, int frame_range)
	{
		return cachehints == CACHE_GET_MTMODE ? MT_MULTI_INSTANCE : 0;
	}
	~vsTBilateral();
};

#define MINS 0.00000000000001

static double kernelValue(double x, double sigma, int kernel)
{
	switch (kernel)
	{
	case AndrewsWave: // Andrews' wave
		if (x <= sigma)
			return ((sin((M_PI * x) / sigma) * sigma) / M_PI);
		return 0.0;
	case ElFallahFord: // El Fallah Ford
		return (1.0 / sqrt(1.0 + ((x * x) / (sigma * sigma))));
	case Gaussian: // Gaussian
		return (exp(-((x * x) / (2.0 * sigma * sigma))));
	case HubersMiniMax: // Huber’s mini-max
		if (x <= sigma)
			return (1.0 / sigma);
		return (1.0 / x);
	case Lorentzian: // Lorentzian
		return (2.0 / (2.0 * sigma * sigma + x * x));
	case TukeyBiWeight: // Tukey bi-weight
		if (x <= sigma)
			return (0.5 * pow((1.0 - ((x * x) / (sigma * sigma))), 2));
		return 0.0;
	case LinearDescent: // Linear descent
		if (x <= sigma)
			return (1.0 - (x / sigma));
		return 0.0;
	case Cosine: // Cosine
		if (x <= sigma)
			return (cos((M_PI * x) / (2.0 * sigma)));
		return 0.0;
	case Flat: // Flat
		if (x <= sigma)
			return (1.0 / sigma);
		return 0.0;
	case Inverse: // Inverse
		if (x <= sigma)
		{
			if (x != 0.0)
				return (1.0 / x);
			return 1.0;
		}
		return 0.0;
	}
	return 0.0;
}

// Singular Value Decomposition routine
// taken from numerical recipes in C.

static double pythag(double a, double b)
{
	double at = fabs(a), bt = fabs(b), ct;
	if (at > bt)
	{
		ct = bt / at;
		return at * sqrt(1.0 + ct * ct);
	}

	if (bt > 0.0)
	{
		ct = at / bt;
		return bt * sqrt(1.0 + ct * ct);
	}

	return 0.0;
}

static void svdcmp(double* a, double* w, double* v)
{
	int flag, i, its, j, jj, k, l, nm;
	double c, f, h, s, x, y, z;
	double anorm = 0.0, g = 0.0, scale = 0.0;
	double rv1[3];
	for (i = 0; i < 3; i++)
	{
		l = i + 1;
		rv1[i] = scale * g;
		g = s = scale = 0.0;

		if (i < 3)
		{
			for (k = i; k < 3; k++)
				scale += fabs(a[k * 3 + i]);

			if (scale)
			{
				for (k = i; k < 3; k++)
				{
					a[k * 3 + i] = a[k * 3 + i] / scale;
					s += a[k * 3 + i] * a[k * 3 + i];
				}

				f = a[i * 3 + i];
				g = -copysign(sqrt(s), f);
				h = f * g - s;
				a[i * 3 + i] = f - g;

				if (i != 2)
				{
					for (j = l; j < 3; j++)
					{
						for (s = 0.0, k = i; k < 3; k++)
							s += a[k * 3 + i] * a[k * 3 + j];
						f = s / h;

						for (k = i; k < 3; k++)
							a[k * 3 + j] += f * a[k * 3 + i];
					}
				}

				for (k = i; k < 3; k++)
					a[k * 3 + i] = a[k * 3 + i] * scale;
			}
		}

		w[i] = scale * g;
		g = s = scale = 0.0;

		if (i < 3 && i != 2)
		{
			for (k = l; k < 3; k++)
				scale += fabs(a[i * 3 + k]);

			if (scale)
			{
				for (k = l; k < 3; k++)
				{
					a[i * 3 + k] = a[i * 3 + k] / scale;
					s += a[i * 3 + k] * a[i * 3 + k];
				}

				f = a[i * 3 + l];
				g = -copysign(sqrt(s), f);
				h = f * g - s;
				a[i * 3 + l] = f - g;

				for (k = l; k < 3; k++)
					rv1[k] = a[i * 3 + k] / h;
				if (i != 2)
				{
					for (j = l; j < 3; j++)
					{
						for (s = 0.0, k = l; k < 3; k++)
							s += (a[j * 3 + k] * a[i * 3 + k]);

						for (k = l; k < 3; k++)
							a[j * 3 + k] += s * rv1[k];
					}
				}

				for (k = l; k < 3; k++)
					a[i * 3 + k] = a[i * 3 + k] * scale;
			}
		}

		anorm = std::max(anorm, (fabs(w[i]) + fabs(rv1[i])));
	}
	for (i = 2; i >= 0; i--)
	{
		if (i < 2)
		{
			if (g)
			{
				for (j = l; j < 3; j++)
					v[j * 3 + i] = a[i * 3 + j] / a[i * 3 + l] / g;

				for (j = l; j < 3; j++)
				{
					for (s = 0.0, k = l; k < 3; k++)
						s += (a[i * 3 + k] * v[k * 3 + j]);

					for (k = l; k < 3; k++)
						v[k * 3 + j] += s * v[k * 3 + i];
				}
			}

			for (j = l; j < 3; j++)
				v[i * 3 + j] = v[j * 3 + i] = 0.0;
		}

		v[i * 3 + i] = 1.0;
		g = rv1[i];
		l = i;
	}

	for (i = 2; i >= 0; i--)
	{
		l = i + 1;
		g = w[i];

		if (i < 2)
			for (j = l; j < 3; j++)
				a[i * 3 + j] = 0.0;

		if (g)
		{
			g = 1.0 / g;

			if (i != 2)
			{
				for (j = l; j < 3; j++)
				{
					for (s = 0.0, k = l; k < 3; k++)
						s += (a[k * 3 + i] * a[k * 3 + j]);
					f = (s / a[i * 3 + i]) * g;

					for (k = i; k < 3; k++)
						a[k * 3 + j] += f * a[k * 3 + i];
				}
			}

			for (j = i; j < 3; j++)
				a[j * 3 + i] = a[j * 3 + i] * g;
		}

		else
		{
			for (j = i; j < 3; j++)
				a[j * 3 + i] = 0.0;
		}
		++a[i * 3 + i];
	}

	for (k = 2; k >= 0; k--)
	{
		for (its = 0; its < 30; its++)
		{
			flag = 1;

			for (l = k; l >= 0; l--)
			{
				nm = l - 1;

				if (fabs(rv1[l]) + anorm == anorm)
				{
					flag = 0;
					break;
				}

				if (fabs(w[nm]) + anorm == anorm)
					break;
			}

			if (flag)
			{
				c = 0.0;
				s = 1.0;

				for (i = l; i <= k; i++)
				{
					f = s * rv1[i];

					if (fabs(f) + anorm != anorm)
					{
						g = w[i];
						h = pythag(f, g);
						w[i] = h;
						h = 1.0 / h;
						c = g * h;
						s = (-f * h);

						for (j = 0; j < 3; j++)
						{
							y = a[j * 3 + nm];
							z = a[j * 3 + i];
							a[j * 3 + nm] = y * c + z * s;
							a[j * 3 + i] = z * c - y * s;
						}
					}
				}
			}

			z = w[k];

			if (l == k)
			{
				if (z < 0.0)
				{
					w[k] = -z;

					for (j = 0; j < 3; j++)
						v[j * 3 + k] = -v[j * 3 + k];
				}
				break;
			}

			x = w[l];
			nm = k - 1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			g = pythag(f, 1.0);
			f = ((x - z) * (x + z) + h * ((y / (f + copysign(g, f))) - h)) / x;
			c = s = 1.0;

			for (j = l; j <= nm; j++)
			{
				i = j + 1;
				g = rv1[i];
				y = w[i];
				h = s * g;
				g = c * g;
				z = pythag(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y = y * c;

				for (jj = 0; jj < 3; jj++)
				{
					x = v[jj * 3 + j];
					z = v[jj * 3 + i];
					v[jj * 3 + j] = x * c + z * s;
					v[jj * 3 + i] = z * c - x * s;
				}

				z = pythag(f, h);
				w[j] = z;

				if (z)
				{
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}

				f = (c * g) + (s * y);
				x = (c * y) - (s * g);

				for (jj = 0; jj < 3; jj++)
				{
					y = a[jj * 3 + j];
					z = a[jj * 3 + i];
					a[jj * 3 + j] = y * c + z * s;
					a[jj * 3 + i] = z * c - y * s;
				}
			}

			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
}

static int mlre(double* yi, double* wi, int lw, int lh, int cx, int cy, int radius, int diameter, int pixel_max)
{
	wi += static_cast<int64_t>(cy) * diameter;
	yi += static_cast<int64_t>(cy) * diameter;

	const int la = lw * lh;
	const int lax2 = la * 2;
	const int la2 = la * la;

	double* xlr = reinterpret_cast<double*>(malloc(la * static_cast<int64_t>(3) * sizeof(double)));
	double* wlr = reinterpret_cast<double*>(malloc(la2 * sizeof(double)));
	double* ylr = reinterpret_cast<double*>(malloc(la * sizeof(double)));
	double* wxlr = reinterpret_cast<double*>(malloc(la * static_cast<int64_t>(3) * sizeof(double)));
	double* xtlr = reinterpret_cast<double*>(malloc(la * static_cast<int64_t>(3) * sizeof(double)));
	double* wylr = reinterpret_cast<double*>(malloc(la * sizeof(double)));
	double xtwx[9], xtwxi[9], xtwy[3], blr[3], wjlr[3], vlr[9];

	// compute w and y matrices
	int d = 0, h = 0;
	memset(wlr, 0, la2 * sizeof(double));
	for (int k = 0; k < lh; ++k)
	{
		const int kt = k * diameter;

		for (int j = cx; j < lw + cx; ++j, ++h, d += la + 1)
		{
			wlr[d] = wi[kt + j];
			ylr[h] = yi[kt + j];
		}
	}

	// compute x and x' matrices
	d = 0;

	for (int j = 0; j < lh; ++j)
	{
		const int jt = j * lw * 3;

		for (int k = 0; k < lw; ++k, ++d)
		{
			xlr[jt + k * 3 + 0] = xtlr[d] = 1;
			xlr[jt + k * 3 + 1] = xtlr[d + la] = j;
			xlr[jt + k * 3 + 2] = xtlr[d + lax2] = k;
		}
	}

	// compute w*x matrix
	for (int j = 0; j < la; ++j)
	{
		const int j3 = j * 3;
		const int jl = j * la;

		for (int k = 0; k < 3; ++k)
		{
			wxlr[j3 + k] = 0.0;

			for (int l = 0; l < la; ++l)
				wxlr[j3 + k] += wlr[jl + l] * xlr[l * 3 + k];
		}
	}

	// compute xt*wx matrix
	for (int j = 0; j < 3; ++j)
	{
		const int j3 = j * 3;
		const int jl = j * la;

		for (int k = 0; k < 3; ++k)
		{
			xtwx[j3 + k] = 0.0;

			for (int l = 0; l < la; ++l)
				xtwx[j3 + k] += xtlr[jl + l] * wxlr[l * 3 + k];
		}
	}

	// compute svd of xtwx = U*WJ*V'
	svdcmp(xtwx, wjlr, vlr);

	// compute wj inverse + zero small wj's
	for (int i = 0; i < 3; ++i)
	{
		if (fabs(wjlr[i]) <= FLT_EPSILON)
			wjlr[i] = 0;
		else
			wjlr[i] = 1.0 / wjlr[i];
	}

	// compute wj^-1 * u'
	for (int j = 0; j < 3; ++j)
	{
		const int j3 = j * 3;

		for (int k = j; k < 3; ++k)
		{
			double temp = xtwx[j3 + k];
			xtwx[j3 + k] = xtwx[k * 3 + j] * wjlr[j];
			xtwx[k * 3 + j] = temp * wjlr[k];
		}
	}

	// compute xtwxi
	for (int j = 0; j < 3; ++j)
	{
		const int j3 = j * 3;

		for (int k = 0; k < 3; ++k)
		{
			xtwxi[j3 + k] = 0.0;

			for (int l = 0; l < 3; ++l)
				xtwxi[j3 + k] += vlr[j * 3 + l] * xtwx[l * 3 + k];
		}
	}

	// compute wy matrix
	for (int j = 0; j < la; ++j)
	{
		const int jl = j * la;
		wylr[j] = 0.0;

		for (int l = 0; l < la; ++l)
			wylr[j] += wlr[jl + l] * ylr[l];
	}

	// compute xtwy matrix
	for (int j = 0; j < 3; ++j)
	{
		const int jl = j * la;
		xtwy[j] = 0.0;

		for (int l = 0; l < la; ++l)
			xtwy[j] += xtlr[jl + l] * wylr[l];
	}

	// compute b matrix
	for (int j = 0; j < 3; ++j)
	{
		const int j3 = j * 3;
		blr[j] = 0.0;

		for (int l = 0; l < 3; ++l)
			blr[j] += xtwxi[j3 + l] * xtwy[l];
	}

	free(xlr);
	free(wlr);
	free(ylr);
	free(wxlr);
	free(xtlr);
	free(wylr);

	return std::min(std::max(int(blr[0] + blr[1] * (radius - static_cast<int64_t>(cy)) + blr[2] * (radius - static_cast<int64_t>(cx)) + 0.5), 0), pixel_max);
}

template <typename PixelType>
void vsTBilateral::ProcessFrameD2_Mean(const uint8_t* srcp_, const uint8_t* tp_, uint8_t* dstp_, int src_pitch, int tp_pitch, int dst_pitch, const int width, const int height, const int plane, IScriptEnvironment* env)
{
	const int pixel_max = vsTBilateral::pixel_max;

	src_pitch /= sizeof(PixelType);
	tp_pitch /= sizeof(PixelType);
	dst_pitch /= sizeof(PixelType);
	const PixelType* srcp = reinterpret_cast<const PixelType*>(srcp_);
	const PixelType* tp = reinterpret_cast<const PixelType*>(tp_);
	PixelType* dstp = reinterpret_cast<PixelType*>(dstp_);

	const int diameter = vsTBilateral::diameter[plane];
	const int radius = diameter >> 1;
	int stopy = radius;
	int startyr = radius * diameter;
	const double* spatialWeights = vsTBilateral::spatialWeights[plane];
	const double* diffWeights = vsTBilateral::diffWeights[plane] + pixel_max * static_cast<int64_t>(2);
	const PixelType* srcp_saved = srcp;
	int starty = 0;
	const int midP = width - radius;
	const int midPY = height - radius;

	int y;
	for (y = 0; y < radius; ++y, startyr -= diameter, ++stopy)
	{
		const PixelType* srcpT_saved = srcp_saved + static_cast<int64_t>(starty) * src_pitch;
		int startx = 0;
		int startxr = radius;
		int stopx = radius;

		int x;
		for (x = 0; x < radius; ++x, --startxr, ++stopx)
		{
			double weightedSum = 0.0;
			double sumOfWeights = 0.0;
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					weightedSum += srcpT[v] * weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = int((weightedSum / sumOfWeights) + 0.5);
			else
				dstp[x] = srcp[x];
		}

		for (; x < midP; ++x, ++startx, ++stopx)
		{
			double weightedSum = 0.0;
			double sumOfWeights = 0.0;
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					weightedSum += srcpT[v] * weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = int((weightedSum / sumOfWeights) + 0.5);
			else
				dstp[x] = srcp[x];
		}

		for (--stopx; x < width; ++x, ++startx)
		{
			double weightedSum = 0.0;
			double sumOfWeights = 0.0;
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					weightedSum += srcpT[v] * weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = int((weightedSum / sumOfWeights) + 0.5);
			else
				dstp[x] = srcp[x];
		}

		srcp += src_pitch;
		dstp += dst_pitch;
		tp += tp_pitch;
	}

	for (; y < midPY; ++y, ++starty, ++stopy)
	{
		const PixelType* srcpT_saved = srcp_saved + static_cast<int64_t>(starty) * src_pitch;
		const PixelType* srcpT2_saved = srcp_saved + static_cast<int64_t>(stopy) * src_pitch;
		int startx = 0;
		int startxr = radius;
		int stopx = radius;

		int x;
		for (x = 0; x < radius; ++x, --startxr, ++stopx)
		{
			double weightedSum = 0.0;
			double sumOfWeights = 0.0;
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					weightedSum += srcpT[v] * weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = int((weightedSum / sumOfWeights) + 0.5);
			else
				dstp[x] = srcp[x];
		}

		for (; x < midP; ++x, ++startx, ++stopx) // free of all boundaries
		{
			double weightedSum = 0.0;
			double sumOfWeights = 0.0;
			const PixelType* srcpT = srcpT_saved;
			const PixelType* srcpT2 = srcpT2_saved;
			const int cP = tp[x] << 1;
			int w = 0;

			for (int u = starty; u <= stopy; ++u)
			{
				int b = stopx;

				for (int v = startx; v <= stopx; ++v, --b, ++w)
				{
					const double weight = spatialWeights[w] * diffWeights[cP - srcpT[v] - srcpT2[b]];
					weightedSum += srcpT[v] * weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
				srcpT2 -= src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = int((weightedSum / sumOfWeights) + 0.5);
			else
				dstp[x] = srcp[x];
		}

		for (--stopx; x < width; ++x, ++startx)
		{
			double weightedSum = 0.0;
			double sumOfWeights = 0.0;
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight =
						spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					weightedSum += srcpT[v] * weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = int((weightedSum / sumOfWeights) + 0.5);
			else
				dstp[x] = srcp[x];
		}

		srcp += src_pitch;
		dstp += dst_pitch;
		tp += tp_pitch;
	}

	for (--stopy; y < height; ++y, ++starty)
	{
		const PixelType* srcpT_saved = srcp_saved + static_cast<int64_t>(starty) * src_pitch;
		int startx = 0;
		int startxr = radius;
		int stopx = radius;

		int x;
		for (x = 0; x < radius; ++x, --startxr, ++stopx)
		{
			double weightedSum = 0.0;
			double sumOfWeights = 0.0;
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight =
						spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					weightedSum += srcpT[v] * weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = int((weightedSum / sumOfWeights) + 0.5);
			else
				dstp[x] = srcp[x];
		}

		for (; x < midP; ++x, ++startx, ++stopx)
		{
			double weightedSum = 0.0;
			double sumOfWeights = 0.0;
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight =
						spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					weightedSum += srcpT[v] * weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = int((weightedSum / sumOfWeights) + 0.5);
			else
				dstp[x] = srcp[x];
		}

		for (--stopx; x < width; ++x, ++startx)
		{
			double weightedSum = 0.0;
			double sumOfWeights = 0.0;
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight =
						spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					weightedSum += srcpT[v] * weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = int((weightedSum / sumOfWeights) + 0.5);
			else
				dstp[x] = srcp[x];
		}

		srcp += src_pitch;
		dstp += dst_pitch;
		tp += tp_pitch;
	}
}

template <typename PixelType>
void vsTBilateral::ProcessFrameD2_MLR(const uint8_t* srcp_, const uint8_t* tp_, uint8_t* dstp_, int src_pitch, int tp_pitch, int dst_pitch, const int width, const int height, const int plane, IScriptEnvironment* env)
{
	const int pixel_max = vsTBilateral::pixel_max;

	src_pitch /= sizeof(PixelType);
	tp_pitch /= sizeof(PixelType);
	dst_pitch /= sizeof(PixelType);
	const PixelType* srcp = reinterpret_cast<const PixelType*>(srcp_);
	const PixelType* tp = reinterpret_cast<const PixelType*>(tp_);
	PixelType* dstp = reinterpret_cast<PixelType*>(dstp_);

	const int diameter = vsTBilateral::diameter[plane];
	const int radius = diameter >> 1;
	int stopy = radius;
	int startyr = radius * diameter;
	const double* spatialWeights = vsTBilateral::spatialWeights[plane];
	const double* diffWeights = vsTBilateral::diffWeights[plane] + pixel_max * static_cast<int64_t>(2);

	const size_t wda = static_cast<int64_t>(diameter) * diameter * sizeof(double);

	double* pixels = reinterpret_cast<double*>(_aligned_malloc(wda, 16));
	double* weights = reinterpret_cast<double*>(_aligned_malloc(wda, 16));

	const PixelType* srcp_saved = srcp;
	int starty = 0;
	const int midP = width - radius;
	const int midPY = height - radius;

	int y;
	for (y = 0; y < radius; ++y, startyr -= diameter, ++stopy)
	{
		const PixelType* srcpT_saved = srcp_saved + static_cast<int64_t>(starty) * src_pitch;
		int startx = 0;
		int startxr = radius;
		int stopx = radius;

		int x;
		for (x = 0; x < radius; ++x, --startxr, ++stopx)
		{
			double sumOfWeights = 0.0;
			memset(pixels, 0, wda);
			memset(weights, 0, wda);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					pixels[w + b] = srcpT[v];
					weights[w + b] = weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = mlre(pixels, weights, stopx - startx + 1,
					stopy - starty + 1, startxr, radius + starty - y, radius, diameter, pixel_max);
			else
				dstp[x] = srcp[x];
		}

		for (; x < midP; ++x, ++startx, ++stopx)
		{
			double sumOfWeights = 0.0;
			memset(pixels, 0, wda);
			memset(weights, 0, wda);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					pixels[w + b] = srcpT[v];
					weights[w + b] = weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = mlre(pixels, weights, diameter,
					stopy - starty + 1, 0, radius + starty - y, radius, diameter, pixel_max);
			else
				dstp[x] = srcp[x];
		}

		for (--stopx; x < width; ++x, ++startx)
		{
			double sumOfWeights = 0.0;
			memset(pixels, 0, wda);
			memset(weights, 0, wda);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					pixels[w + b] = srcpT[v];
					weights[w + b] = weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = mlre(pixels, weights, stopx - startx + 1,
					stopy - starty + 1, startxr, radius + starty - y, radius, diameter, pixel_max);
			else
				dstp[x] = srcp[x];
		}

		srcp += src_pitch;
		dstp += dst_pitch;
		tp += tp_pitch;
	}

	for (; y < midPY; ++y, ++starty, ++stopy)
	{
		const PixelType* srcpT_saved = srcp_saved + static_cast<int64_t>(starty) * src_pitch;
		const PixelType* srcpT2_saved = srcp_saved + static_cast<int64_t>(stopy) * src_pitch;
		int startx = 0;
		int startxr = radius;
		int stopx = radius;

		int x;
		for (x = 0; x < radius; ++x, --startxr, ++stopx)
		{
			double sumOfWeights = 0.0;
			memset(pixels, 0, wda);
			memset(weights, 0, wda);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					pixels[w + b] = srcpT[v];
					weights[w + b] = weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = mlre(pixels, weights, stopx - startx + 1,
					diameter, startxr, 0, radius, diameter, pixel_max);
			else
				dstp[x] = srcp[x];
		}

		for (; x < midP; ++x, ++startx, ++stopx) // free of all boundaries
		{
			double sumOfWeights = 0.0;
			memset(pixels, 0, wda);
			memset(weights, 0, wda);
			const PixelType* srcpT = srcpT_saved;
			const PixelType* srcpT2 = srcpT2_saved;
			const int cP = tp[x] << 1;
			int w = 0;

			for (int u = starty; u <= stopy; ++u)
			{
				int b = stopx;

				for (int v = startx; v <= stopx; ++v, --b, ++w)
				{
					const double weight = spatialWeights[w] * diffWeights[cP - srcpT[v] - srcpT2[b]];
					pixels[w] = srcpT[v];
					weights[w] = weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
				srcpT2 -= src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = mlre(pixels, weights, diameter,
					diameter, 0, 0, radius, diameter, pixel_max);
			else
				dstp[x] = srcp[x];
		}

		for (--stopx; x < width; ++x, ++startx)
		{
			double sumOfWeights = 0.0;
			memset(pixels, 0, wda);
			memset(weights, 0, wda);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					pixels[w + b] = srcpT[v];
					weights[w + b] = weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = mlre(pixels, weights, stopx - startx + 1,
					diameter, startxr, 0, radius, diameter, pixel_max);
			else
				dstp[x] = srcp[x];
		}

		srcp += src_pitch;
		dstp += dst_pitch;
		tp += tp_pitch;
	}

	for (--stopy; y < height; ++y, ++starty)
	{
		const PixelType* srcpT_saved = srcp_saved + static_cast<int64_t>(starty) * src_pitch;
		int startx = 0;
		int startxr = radius;
		int stopx = radius;

		int x;
		for (x = 0; x < radius; ++x, --startxr, ++stopx)
		{
			double sumOfWeights = 0.0;
			memset(pixels, 0, wda);
			memset(weights, 0, wda);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					pixels[w + b] = srcpT[v];
					weights[w + b] = weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = mlre(pixels, weights, stopx - startx + 1,
					stopy - starty + 1, startxr, radius + starty - y, radius, diameter, pixel_max);
			else
				dstp[x] = srcp[x];
		}

		for (; x < midP; ++x, ++startx, ++stopx)
		{
			double sumOfWeights = 0.0;
			memset(pixels, 0, wda);
			memset(weights, 0, wda);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					pixels[w + b] = srcpT[v];
					weights[w + b] = weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = mlre(pixels, weights, diameter,
					stopy - starty + 1, 0, radius + starty - y, radius, diameter, pixel_max);
			else
				dstp[x] = srcp[x];
		}

		for (--stopx; x < width; ++x, ++startx)
		{
			double sumOfWeights = 0.0;
			memset(pixels, 0, wda);
			memset(weights, 0, wda);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					pixels[w + b] = srcpT[v];
					weights[w + b] = weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = mlre(pixels, weights, stopx - startx + 1,
					stopy - starty + 1, startxr, radius + starty - y, radius, diameter, pixel_max);
			else
				dstp[x] = srcp[x];
		}

		srcp += src_pitch;
		dstp += dst_pitch;
		tp += tp_pitch;
	}

	_aligned_free(pixels);
	_aligned_free(weights);
}

template <typename PixelType>
void vsTBilateral::ProcessFrameD2_Med(const uint8_t* srcp_, const uint8_t* tp_, uint8_t* dstp_, int src_pitch, int tp_pitch, int dst_pitch, const int width, const int height, const int plane, IScriptEnvironment* env)
{
	const bool cw = resType_ == CWMedian;

	const int pixel_max = vsTBilateral::pixel_max;

	size_t medAsize = (pixel_max + static_cast<int64_t>(1)) * sizeof(double);
	double* medA = reinterpret_cast<double*>(_aligned_malloc(medAsize, 16));

	src_pitch /= sizeof(PixelType);
	tp_pitch /= sizeof(PixelType);
	dst_pitch /= sizeof(PixelType);
	const PixelType* srcp = reinterpret_cast<const PixelType*>(srcp_);
	const PixelType* tp = reinterpret_cast<const PixelType*>(tp_);
	PixelType* dstp = reinterpret_cast<PixelType*>(dstp_);

	const int diameter = vsTBilateral::diameter[plane];
	const int radius = diameter >> 1;
	int stopy = radius;
	int startyr = radius * diameter;
	const double* spatialWeights = vsTBilateral::spatialWeights[plane];
	const double* diffWeights = vsTBilateral::diffWeights[plane] + pixel_max * static_cast<int64_t>(2);
	const PixelType* srcp_saved = srcp;
	int starty = 0;
	const int midP = width - radius;
	const int midPY = height - radius;
	const int mid = diameter * radius + radius;
	const double cw_weight = spatialWeights[mid] * diffWeights[-pixel_max * 2] * (diameter - static_cast<int64_t>(1));
	int y;

	for (y = 0; y < radius; ++y, startyr -= diameter, ++stopy)
	{
		const PixelType* srcpT_saved = srcp_saved + static_cast<int64_t>(starty) * src_pitch;
		int startx = 0;
		int startxr = radius;
		int stopx = radius;

		int x;
		for (x = 0; x < radius; ++x, --startxr, ++stopx)
		{
			double sumOfWeights = 0.0;
			double sum = 0.0;
			memset(medA, 0, medAsize);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					medA[srcpT[v]] += weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
			{
				if (cw)
				{
					medA[tp[x]] += cw_weight;
					sumOfWeights += cw_weight;
				}

				sumOfWeights *= 0.5;
				int ws = 0;

				while (sum <= sumOfWeights)
				{
					sum += medA[ws];
					++ws;
				}

				dstp[x] = ws - 1;
			}
			else
				dstp[x] = srcp[x];
		}

		for (; x < midP; ++x, ++startx, ++stopx)
		{
			double sumOfWeights = 0.0;
			double sum = 0.0;
			memset(medA, 0, medAsize);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					medA[srcpT[v]] += weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
			{
				if (cw)
				{
					medA[tp[x]] += cw_weight;
					sumOfWeights += cw_weight;
				}

				sumOfWeights *= 0.5;
				int ws = 0;

				while (sum <= sumOfWeights)
				{
					sum += medA[ws];
					++ws;
				}

				dstp[x] = ws - 1;
			}
			else
				dstp[x] = srcp[x];
		}

		for (--stopx; x < width; ++x, ++startx)
		{
			double sumOfWeights = 0.0;
			double sum = 0.0;
			memset(medA, 0, medAsize);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					medA[srcpT[v]] += weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
			{
				if (cw)
				{
					medA[tp[x]] += cw_weight;
					sumOfWeights += cw_weight;
				}

				sumOfWeights *= 0.5;
				int ws = 0;

				while (sum <= sumOfWeights)
				{
					sum += medA[ws];
					++ws;
				}

				dstp[x] = ws - 1;
			}
			else
				dstp[x] = srcp[x];
		}

		srcp += src_pitch;
		dstp += dst_pitch;
		tp += tp_pitch;
	}

	for (; y < midPY; ++y, ++starty, ++stopy)
	{
		const PixelType* srcpT_saved = srcp_saved + static_cast<int64_t>(starty) * src_pitch;
		const PixelType* srcpT2_saved = srcp_saved + static_cast<int64_t>(stopy) * src_pitch;
		int startx = 0;
		int startxr = radius;
		int stopx = radius;

		int x;
		for (x = 0; x < radius; ++x, --startxr, ++stopx)
		{
			double sumOfWeights = 0.0;
			double sum = 0.0;
			memset(medA, 0, medAsize);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					medA[srcpT[v]] += weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
			{
				if (cw)
				{
					medA[tp[x]] += cw_weight;
					sumOfWeights += cw_weight;
				}

				sumOfWeights *= 0.5;
				int ws = 0;

				while (sum <= sumOfWeights)
				{
					sum += medA[ws];
					++ws;
				}

				dstp[x] = ws - 1;
			}
			else
				dstp[x] = srcp[x];
		}

		for (; x < midP; ++x, ++startx, ++stopx) // free of all boundaries
		{
			double sumOfWeights = 0.0;
			double sum = 0.0;
			memset(medA, 0, medAsize);
			const PixelType* srcpT = srcpT_saved;
			const PixelType* srcpT2 = srcpT2_saved;
			const int cP = tp[x] << 1;
			int w = 0;

			for (int u = starty; u <= stopy; ++u)
			{
				int b = stopx;

				for (int v = startx; v <= stopx; ++v, --b, ++w)
				{
					const double weight = spatialWeights[w] * diffWeights[cP - srcpT[v] - srcpT2[b]];
					medA[srcpT[v]] += weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
				srcpT2 -= src_pitch;
			}

			if (sumOfWeights >= MINS)
			{
				if (cw)
				{
					medA[tp[x]] += cw_weight;
					sumOfWeights += cw_weight;
				}

				sumOfWeights *= 0.5;
				int ws = 0;

				while (sum <= sumOfWeights)
				{
					sum += medA[ws];
					++ws;
				}

				dstp[x] = ws - 1;
			}
			else
				dstp[x] = srcp[x];
		}

		for (--stopx; x < width; ++x, ++startx)
		{
			double sumOfWeights = 0.0;
			double sum = 0.0;
			memset(medA, 0, medAsize);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					medA[srcpT[v]] += weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
			{
				if (cw)
				{
					medA[tp[x]] += cw_weight;
					sumOfWeights += cw_weight;
				}

				sumOfWeights *= 0.5;
				int ws = 0;

				while (sum <= sumOfWeights)
				{
					sum += medA[ws];
					++ws;
				}

				dstp[x] = ws - 1;
			}
			else
				dstp[x] = srcp[x];
		}

		srcp += src_pitch;
		dstp += dst_pitch;
		tp += tp_pitch;
	}

	for (--stopy; y < height; ++y, ++starty)
	{
		const PixelType* srcpT_saved = srcp_saved + static_cast<int64_t>(starty) * src_pitch;
		int startx = 0;
		int startxr = radius;
		int stopx = radius;

		int x;
		for (x = 0; x < radius; ++x, --startxr, ++stopx)
		{
			double sumOfWeights = 0.0;
			double sum = 0.0;
			memset(medA, 0, medAsize);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					medA[srcpT[v]] += weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
			{
				if (cw)
				{
					medA[tp[x]] += cw_weight;
					sumOfWeights += cw_weight;
				}

				sumOfWeights *= 0.5;
				int ws = 0;

				while (sum <= sumOfWeights)
				{
					sum += medA[ws];
					++ws;
				}

				dstp[x] = ws - 1;
			}
			else
				dstp[x] = srcp[x];
		}

		for (; x < midP; ++x, ++startx, ++stopx)
		{
			double sumOfWeights = 0.0;
			double sum = 0.0;
			memset(medA, 0, medAsize);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					medA[srcpT[v]] += weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
			{
				if (cw)
				{
					medA[tp[x]] += cw_weight;
					sumOfWeights += cw_weight;
				}

				sumOfWeights *= 0.5;
				int ws = 0;

				while (sum <= sumOfWeights)
				{
					sum += medA[ws];
					++ws;
				}

				dstp[x] = ws - 1;
			}
			else
				dstp[x] = srcp[x];
		}

		for (--stopx; x < width; ++x, ++startx)
		{
			double sumOfWeights = 0.0;
			double sum = 0.0;
			memset(medA, 0, medAsize);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;
				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[(cP - srcpT[v]) << 1];
					medA[srcpT[v]] += weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
			{
				if (cw)
				{
					medA[tp[x]] += cw_weight;
					sumOfWeights += cw_weight;
				}

				sumOfWeights *= 0.5;
				int ws = 0;

				while (sum <= sumOfWeights)
				{
					sum += medA[ws];
					++ws;
				}

				dstp[x] = ws - 1;
			}
			else
				dstp[x] = srcp[x];
		}

		srcp += src_pitch;
		dstp += dst_pitch;
		tp += tp_pitch;
	}

	_aligned_free(medA);
}

template <typename PixelType>
void vsTBilateral::ProcessFrameD1_Mean(const uint8_t* srcp_, const uint8_t* tp_, uint8_t* dstp_, int src_pitch, int tp_pitch, int dst_pitch, const int width, const int height, const int plane, IScriptEnvironment* env)
{
	const int pixel_max = vsTBilateral::pixel_max;

	src_pitch /= sizeof(PixelType);
	tp_pitch /= sizeof(PixelType);
	dst_pitch /= sizeof(PixelType);
	const PixelType* srcp = reinterpret_cast<const PixelType*>(srcp_);
	const PixelType* tp = reinterpret_cast<const PixelType*>(tp_);
	PixelType* dstp = reinterpret_cast<PixelType*>(dstp_);

	const int diameter = vsTBilateral::diameter[plane];
	const int radius = diameter >> 1;
	int stopy = radius;
	int startyr = radius * diameter;
	const double* spatialWeights = vsTBilateral::spatialWeights[plane];
	const double* diffWeights = vsTBilateral::diffWeights[plane] + pixel_max;
	const PixelType* srcp_saved = srcp;
	int starty = 0;
	const int midP = width - radius;
	const int midPY = height - radius;

	int y;
	for (y = 0; y < radius; ++y, startyr -= diameter, ++stopy)
	{
		const PixelType* srcpT_saved = srcp_saved + static_cast<int64_t>(starty) * src_pitch;
		int startx = 0;
		int startxr = radius;
		int stopx = radius;

		int x;
		for (x = 0; x < radius; ++x, --startxr, ++stopx)
		{
			double weightedSum = 0.0;
			double sumOfWeights = 0.0;
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					weightedSum += srcpT[v] * weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = int((weightedSum / sumOfWeights) + 0.5);
			else
				dstp[x] = srcp[x];
		}

		for (; x < midP; ++x, ++startx, ++stopx)
		{
			double weightedSum = 0.0;
			double sumOfWeights = 0.0;
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					weightedSum += srcpT[v] * weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = int((weightedSum / sumOfWeights) + 0.5);
			else
				dstp[x] = srcp[x];
		}

		for (--stopx; x < width; ++x, ++startx)
		{
			double weightedSum = 0.0;
			double sumOfWeights = 0.0;
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					weightedSum += srcpT[v] * weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = int((weightedSum / sumOfWeights) + 0.5);
			else
				dstp[x] = srcp[x];
		}

		srcp += src_pitch;
		dstp += dst_pitch;
		tp += tp_pitch;
	}

	for (; y < midPY; ++y, ++starty, ++stopy)
	{
		const PixelType* srcpT_saved = srcp_saved + static_cast<int64_t>(starty) * src_pitch;
		int startx = 0;
		int startxr = radius;
		int stopx = radius;

		int x;
		for (x = 0; x < radius; ++x, --startxr, ++stopx)
		{
			double weightedSum = 0.0;
			double sumOfWeights = 0.0;
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					weightedSum += srcpT[v] * weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = int((weightedSum / sumOfWeights) + 0.5);
			else
				dstp[x] = srcp[x];
		}

		for (; x < midP; ++x, ++startx, ++stopx) // free of all boundaries
		{
			double weightedSum = 0.0;
			double sumOfWeights = 0.0;
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = 0;

			for (int u = starty; u <= stopy; ++u)
			{
				for (int v = startx; v <= stopx; ++v, ++w)
				{
					const double weight = spatialWeights[w] * diffWeights[cP - srcpT[v]];
					weightedSum += srcpT[v] * weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = int((weightedSum / sumOfWeights) + 0.5);
			else
				dstp[x] = srcp[x];
		}

		for (--stopx; x < width; ++x, ++startx)
		{
			double weightedSum = 0.0;
			double sumOfWeights = 0.0;
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					weightedSum += srcpT[v] * weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = int((weightedSum / sumOfWeights) + 0.5);
			else
				dstp[x] = srcp[x];
		}

		srcp += src_pitch;
		dstp += dst_pitch;
		tp += tp_pitch;
	}

	for (--stopy; y < height; ++y, ++starty)
	{
		const PixelType* srcpT_saved = srcp_saved + static_cast<int64_t>(starty) * src_pitch;
		int startx = 0;
		int startxr = radius;
		int stopx = radius;

		int x;
		for (x = 0; x < radius; ++x, --startxr, ++stopx)
		{
			double weightedSum = 0.0;
			double sumOfWeights = 0.0;
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					weightedSum += srcpT[v] * weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = int((weightedSum / sumOfWeights) + 0.5);
			else
				dstp[x] = srcp[x];
		}

		for (; x < midP; ++x, ++startx, ++stopx)
		{
			double weightedSum = 0.0;
			double sumOfWeights = 0.0;
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					weightedSum += srcpT[v] * weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = int((weightedSum / sumOfWeights) + 0.5);
			else
				dstp[x] = srcp[x];
		}

		for (--stopx; x < width; ++x, ++startx)
		{
			double weightedSum = 0.0;
			double sumOfWeights = 0.0;
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					weightedSum += srcpT[v] * weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = int((weightedSum / sumOfWeights) + 0.5);
			else
				dstp[x] = srcp[x];
		}

		srcp += src_pitch;
		dstp += dst_pitch;
		tp += tp_pitch;
	}
}

template <typename PixelType>
void vsTBilateral::ProcessFrameD1_MLR(const uint8_t* srcp_, const uint8_t* tp_, uint8_t* dstp_, int src_pitch, int tp_pitch, int dst_pitch, const int width, const int height, const int plane, IScriptEnvironment* env)
{
	const int pixel_max = vsTBilateral::pixel_max;

	src_pitch /= sizeof(PixelType);
	tp_pitch /= sizeof(PixelType);
	dst_pitch /= sizeof(PixelType);
	const PixelType* srcp = reinterpret_cast<const PixelType*>(srcp_);
	const PixelType* tp = reinterpret_cast<const PixelType*>(tp_);
	PixelType* dstp = reinterpret_cast<PixelType*>(dstp_);

	const int diameter = vsTBilateral::diameter[plane];
	const int radius = diameter >> 1;
	int stopy = radius;
	int startyr = radius * diameter;
	const double* spatialWeights = vsTBilateral::spatialWeights[plane];
	const double* diffWeights = vsTBilateral::diffWeights[plane] + pixel_max;

	const size_t wda = static_cast<int64_t>(diameter) * diameter * sizeof(double);

	double* pixels = reinterpret_cast<double*>(_aligned_malloc(wda, 16));
	double* weights = reinterpret_cast<double*>(_aligned_malloc(wda, 16));

	const PixelType* srcp_saved = srcp;
	int starty = 0;
	const int midP = width - radius;
	const int midPY = height - radius;

	int y;
	for (y = 0; y < radius; ++y, startyr -= diameter, ++stopy)
	{
		const PixelType* srcpT_saved = srcp_saved + static_cast<int64_t>(starty) * src_pitch;
		int startx = 0;
		int startxr = radius;
		int stopx = radius;

		int x;
		for (x = 0; x < radius; ++x, --startxr, ++stopx)
		{
			double sumOfWeights = 0.0;
			memset(pixels, 0, wda);
			memset(weights, 0, wda);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					pixels[w + b] = srcpT[v];
					weights[w + b] = weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = mlre(pixels, weights, stopx - startx + 1,
					stopy - starty + 1, startxr, radius + starty - y, radius, diameter, pixel_max);
			else
				dstp[x] = srcp[x];
		}

		for (; x < midP; ++x, ++startx, ++stopx)
		{
			double sumOfWeights = 0.0;
			memset(pixels, 0, wda);
			memset(weights, 0, wda);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					pixels[w + b] = srcpT[v];
					weights[w + b] = weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = mlre(pixels, weights, diameter,
					stopy - starty + 1, 0, radius + starty - y, radius, diameter, pixel_max);
			else
				dstp[x] = srcp[x];
		}

		for (--stopx; x < width; ++x, ++startx)
		{
			double sumOfWeights = 0.0;
			memset(pixels, 0, wda);
			memset(weights, 0, wda);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					pixels[w + b] = srcpT[v];
					weights[w + b] = weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = mlre(pixels, weights, stopx - startx + 1,
					stopy - starty + 1, startxr, radius + starty - y, radius, diameter, pixel_max);
			else
				dstp[x] = srcp[x];
		}

		srcp += src_pitch;
		dstp += dst_pitch;
		tp += tp_pitch;
	}

	for (; y < midPY; ++y, ++starty, ++stopy)
	{
		const PixelType* srcpT_saved = srcp_saved + static_cast<int64_t>(starty) * src_pitch;
		int startx = 0;
		int startxr = radius;
		int stopx = radius;

		int x;
		for (x = 0; x < radius; ++x, --startxr, ++stopx)
		{
			double sumOfWeights = 0.0;
			memset(pixels, 0, wda);
			memset(weights, 0, wda);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					pixels[w + b] = srcpT[v];
					weights[w + b] = weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = mlre(pixels, weights, stopx - startx + 1,
					diameter, startxr, 0, radius, diameter, pixel_max);
			else
				dstp[x] = srcp[x];
		}

		for (; x < midP; ++x, ++startx, ++stopx) // free of all boundaries
		{
			double sumOfWeights = 0.0;
			memset(pixels, 0, wda);
			memset(weights, 0, wda);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = 0;

			for (int u = starty; u <= stopy; ++u)
			{
				for (int v = startx; v <= stopx; ++v, ++w)
				{
					const double weight = spatialWeights[w] * diffWeights[cP - srcpT[v]];
					pixels[w] = srcpT[v];
					weights[w] = weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = mlre(pixels, weights, diameter,
					diameter, 0, 0, radius, diameter, pixel_max);
			else
				dstp[x] = srcp[x];
		}

		for (--stopx; x < width; ++x, ++startx)
		{
			double sumOfWeights = 0.0;
			memset(pixels, 0, wda);
			memset(weights, 0, wda);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;
				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					pixels[w + b] = srcpT[v];
					weights[w + b] = weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = mlre(pixels, weights, stopx - startx + 1,
					diameter, startxr, 0, radius, diameter, pixel_max);
			else
				dstp[x] = srcp[x];
		}

		srcp += src_pitch;
		dstp += dst_pitch;
		tp += tp_pitch;
	}

	for (--stopy; y < height; ++y, ++starty)
	{
		const PixelType* srcpT_saved = srcp_saved + static_cast<int64_t>(starty) * src_pitch;
		int startx = 0;
		int startxr = radius;
		int stopx = radius;

		int x;
		for (x = 0; x < radius; ++x, --startxr, ++stopx)
		{
			double sumOfWeights = 0.0;
			memset(pixels, 0, wda);
			memset(weights, 0, wda);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					pixels[w + b] = srcpT[v];
					weights[w + b] = weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = mlre(pixels, weights, stopx - startx + 1,
					stopy - starty + 1, startxr, radius + starty - y, radius, diameter, pixel_max);
			else
				dstp[x] = srcp[x];
		}

		for (; x < midP; ++x, ++startx, ++stopx)
		{
			double sumOfWeights = 0.0;
			memset(pixels, 0, wda);
			memset(weights, 0, wda);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					pixels[w + b] = srcpT[v];
					weights[w + b] = weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = mlre(pixels, weights, diameter,
					stopy - starty + 1, 0, radius + starty - y, radius, diameter, pixel_max);
			else
				dstp[x] = srcp[x];
		}

		for (--stopx; x < width; ++x, ++startx)
		{
			double sumOfWeights = 0.0;
			memset(pixels, 0, wda);
			memset(weights, 0, wda);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					pixels[w + b] = srcpT[v];
					weights[w + b] = weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
				dstp[x] = mlre(pixels, weights, stopx - startx + 1,
					stopy - starty + 1, startxr, radius + starty - y, radius, diameter, pixel_max);
			else
				dstp[x] = srcp[x];
		}

		srcp += src_pitch;
		dstp += dst_pitch;
		tp += tp_pitch;
	}

	_aligned_free(pixels);
	_aligned_free(weights);
}

template <typename PixelType>
void vsTBilateral::ProcessFrameD1_Med(const uint8_t* srcp_, const uint8_t* tp_, uint8_t* dstp_, int src_pitch, int tp_pitch, int dst_pitch, const int width, const int height, const int plane, IScriptEnvironment* env)
{
	const bool cw = resType_ == CWMedian;

	const int pixel_max = vsTBilateral::pixel_max;

	size_t medAsize = (pixel_max + static_cast<int64_t>(1)) * sizeof(double);
	double* medA = reinterpret_cast<double*>(_aligned_malloc(medAsize, 16));

	src_pitch /= sizeof(PixelType);
	tp_pitch /= sizeof(PixelType);
	dst_pitch /= sizeof(PixelType);
	const PixelType* srcp = reinterpret_cast<const PixelType*>(srcp_);
	const PixelType* tp = reinterpret_cast<const PixelType*>(tp_);
	PixelType* dstp = reinterpret_cast<PixelType*>(dstp_);

	const int diameter = vsTBilateral::diameter[plane];
	const int radius = diameter >> 1;
	int stopy = radius;
	int startyr = radius * diameter;
	const double* spatialWeights = vsTBilateral::spatialWeights[plane];
	const double* diffWeights = vsTBilateral::diffWeights[plane] + pixel_max;
	const PixelType* srcp_saved = srcp;
	int starty = 0;
	const int midP = width - radius;
	const int midPY = height - radius;
	const int mid = diameter * radius + radius;
	const double cw_weight = spatialWeights[mid] * diffWeights[-pixel_max] * (diameter - static_cast<int64_t>(1));

	int y;
	for (y = 0; y < radius; ++y, startyr -= diameter, ++stopy)
	{
		const PixelType* srcpT_saved = srcp_saved + static_cast<int64_t>(starty) * src_pitch;
		int startx = 0;
		int startxr = radius;
		int stopx = radius;

		int x;
		for (x = 0; x < radius; ++x, --startxr, ++stopx)
		{
			double sumOfWeights = 0.0;
			double sum = 0.0;
			memset(medA, 0, medAsize);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					medA[srcpT[v]] += weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
			{
				if (cw)
				{
					medA[tp[x]] += cw_weight;
					sumOfWeights += cw_weight;
				}

				sumOfWeights *= 0.5;
				int ws = 0;

				while (sum <= sumOfWeights)
				{
					sum += medA[ws];
					++ws;
				}

				dstp[x] = ws - 1;
			}
			else
				dstp[x] = srcp[x];
		}

		for (; x < midP; ++x, ++startx, ++stopx)
		{
			double sumOfWeights = 0.0;
			double sum = 0.0;
			memset(medA, 0, medAsize);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					medA[srcpT[v]] += weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
			{
				if (cw)
				{
					medA[tp[x]] += cw_weight;
					sumOfWeights += cw_weight;
				}

				sumOfWeights *= 0.5;
				int ws = 0;

				while (sum <= sumOfWeights)
				{
					sum += medA[ws];
					++ws;
				}

				dstp[x] = ws - 1;
			}
			else
				dstp[x] = srcp[x];
		}

		for (--stopx; x < width; ++x, ++startx)
		{
			double sumOfWeights = 0.0;
			double sum = 0.0;
			memset(medA, 0, medAsize);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					medA[srcpT[v]] += weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
			{
				if (cw)
				{
					medA[tp[x]] += cw_weight;
					sumOfWeights += cw_weight;
				}

				sumOfWeights *= 0.5;
				int ws = 0;

				while (sum <= sumOfWeights)
				{
					sum += medA[ws];
					++ws;
				}

				dstp[x] = ws - 1;
			}
			else
				dstp[x] = srcp[x];
		}

		srcp += src_pitch;
		dstp += dst_pitch;
		tp += tp_pitch;
	}

	for (; y < midPY; ++y, ++starty, ++stopy)
	{
		const PixelType* srcpT_saved = srcp_saved + static_cast<int64_t>(starty) * src_pitch;
		int startx = 0;
		int startxr = radius;
		int stopx = radius;

		int x;
		for (x = 0; x < radius; ++x, --startxr, ++stopx)
		{
			double sumOfWeights = 0.0;
			double sum = 0.0;
			memset(medA, 0, medAsize);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					medA[srcpT[v]] += weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
			{
				if (cw)
				{
					medA[tp[x]] += cw_weight;
					sumOfWeights += cw_weight;
				}

				sumOfWeights *= 0.5;
				int ws = 0;

				while (sum <= sumOfWeights)
				{
					sum += medA[ws];
					++ws;
				}

				dstp[x] = ws - 1;
			}
			else
				dstp[x] = srcp[x];
		}

		for (; x < midP; ++x, ++startx, ++stopx) // free of all boundaries
		{
			double sumOfWeights = 0.0;
			double sum = 0.0;
			memset(medA, 0, medAsize);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = 0;

			for (int u = starty; u <= stopy; ++u)
			{
				for (int v = startx; v <= stopx; ++v, ++w)
				{
					const double weight = spatialWeights[w] * diffWeights[cP - srcpT[v]];
					medA[srcpT[v]] += weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
			{
				if (cw)
				{
					medA[tp[x]] += cw_weight;
					sumOfWeights += cw_weight;
				}

				sumOfWeights *= 0.5;
				int ws = 0;

				while (sum <= sumOfWeights)
				{
					sum += medA[ws];
					++ws;
				}

				dstp[x] = ws - 1;
			}
			else
				dstp[x] = srcp[x];
		}

		for (--stopx; x < width; ++x, ++startx)
		{
			double sumOfWeights = 0.0;
			double sum = 0.0;
			memset(medA, 0, medAsize);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					medA[srcpT[v]] += weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
			{
				if (cw)
				{
					medA[tp[x]] += cw_weight;
					sumOfWeights += cw_weight;
				}

				sumOfWeights *= 0.5;
				int ws = 0;

				while (sum <= sumOfWeights)
				{
					sum += medA[ws];
					++ws;
				}

				dstp[x] = ws - 1;
			}
			else
				dstp[x] = srcp[x];
		}

		srcp += src_pitch;
		dstp += dst_pitch;
		tp += tp_pitch;
	}

	for (--stopy; y < height; ++y, ++starty)
	{
		const PixelType* srcpT_saved = srcp_saved + static_cast<int64_t>(starty) * src_pitch;
		int startx = 0;
		int startxr = radius;
		int stopx = radius;

		int x;
		for (x = 0; x < radius; ++x, --startxr, ++stopx)
		{
			double sumOfWeights = 0.0;
			double sum = 0.0;
			memset(medA, 0, medAsize);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					medA[srcpT[v]] += weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
			{
				if (cw)
				{
					medA[tp[x]] += cw_weight;
					sumOfWeights += cw_weight;
				}

				sumOfWeights *= 0.5;
				int ws = 0;

				while (sum <= sumOfWeights)
				{
					sum += medA[ws];
					++ws;
				}

				dstp[x] = ws - 1;
			}
			else
				dstp[x] = srcp[x];
		}

		for (; x < midP; ++x, ++startx, ++stopx)
		{
			double sumOfWeights = 0.0;
			double sum = 0.0;
			memset(medA, 0, medAsize);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					medA[srcpT[v]] += weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
			{
				if (cw)
				{
					medA[tp[x]] += cw_weight;
					sumOfWeights += cw_weight;
				}

				sumOfWeights *= 0.5;
				int ws = 0;

				while (sum <= sumOfWeights)
				{
					sum += medA[ws];
					++ws;
				}

				dstp[x] = ws - 1;
			}
			else
				dstp[x] = srcp[x];
		}

		for (--stopx; x < width; ++x, ++startx)
		{
			double sumOfWeights = 0.0;
			double sum = 0.0;
			memset(medA, 0, medAsize);
			const PixelType* srcpT = srcpT_saved;
			const int cP = tp[x];
			int w = startyr;

			for (int u = starty; u <= stopy; ++u, w += diameter)
			{
				int b = startxr;

				for (int v = startx; v <= stopx; ++v, ++b)
				{
					const double weight = spatialWeights[w + b] * diffWeights[cP - srcpT[v]];
					medA[srcpT[v]] += weight;
					sumOfWeights += weight;
				}

				srcpT += src_pitch;
			}

			if (sumOfWeights >= MINS)
			{
				if (cw)
				{
					medA[tp[x]] += cw_weight;
					sumOfWeights += cw_weight;
				}

				sumOfWeights *= 0.5;
				int ws = 0;

				while (sum <= sumOfWeights)
				{
					sum += medA[ws];
					++ws;
				}

				dstp[x] = ws - 1;
			}
			else
				dstp[x] = srcp[x];
		}

		srcp += src_pitch;
		dstp += dst_pitch;
		tp += tp_pitch;
	}

	_aligned_free(medA);
}

vsTBilateral::vsTBilateral(PClip _child, PClip ppclip, int diameterY, int diameterU, int diameterV, double sdevY, double sdevU, double sdevV, double idevY, double idevU,
	double idevV, double csY, double csU, double csV, bool d2, int kerns, int kerni, int restype, int y, int u, int v, IScriptEnvironment* env)
	: GenericVideoFilter(_child), ppclip_(ppclip), d2_(d2), kerns_(kerns), kerni_(kerni), resType_(restype), y_(y), u_(u), v_(v)
{
	if (kerns_ < AndrewsWave || kerns_ > Inverse)
		env->ThrowError("vsTBilateral: kerns must be between 0 and 9 (inclusive).");
	if (kerni_ < AndrewsWave || kerni_ > Inverse)
		env->ThrowError("vsTBilateral: kerni must be between 0 and 9 (inclusive).");
	if (resType_ < Mean || resType_ > MultipleLinearRegression)
		env->ThrowError("vsTBilateral: restype must be between 0 and 3 (inclusive).");
	if (!vi.IsPlanar() || (vi.IsPlanar() && vi.BitsPerComponent() == 32))
		env->ThrowError("vsTBilateral: clip must be 8..16 bit in planar format.");
	if (vi.width % 2 == 1 || vi.height % 2 == 1)
		env->ThrowError("vsTBilateral: clip's dimensions must be multiples of 2.");
	if (ppclip_)
	{
		const VideoInfo& vi1 = ppclip_->GetVideoInfo();
		if (!vi.IsSameColorspace(vi1) || vi.width != vi1.width || vi.height != vi1.height)
			env->ThrowError("vsTBilateral: ppclip must have the same dimension as the main clip and be the same format.");
		if (vi.num_frames != vi1.num_frames)
			env->ThrowError("vsTBilateral: ppclip's number of frames doesn't match.");
	}

	pixel_max = (1 << vi.BitsPerComponent()) - 1;

	const int planecount = std::min(vi.NumComponents(), 3);
	for (int i = 0; i < planecount; ++i)
	{
		int width[3], height[3];
		switch (i)
		{
		case 0:
		{
			diameter[i] = diameterY;
			sDev[i] = sdevY;
			iDev[i] = idevY;
			cs[i] = csY;
			width[i] = vi.width;
			height[i] = vi.height;
			
			if (vi.IsRGB())
			{
				if (diameter[i] > width[i] || diameter[i] > height[i])
					env->ThrowError("vsTBilateral: diameter must be less than the dimensions of the plane.");

				process[i] = 3;
			}
			else
			{
				if (diameter[i] > width[i] || diameter[i] > height[i])
					env->ThrowError("vsTBilateral: diameterY must be less than the dimensions of plane Y.");

				switch (y_)
				{
				case 3:
					process[i] = 3;
					break;
				case 2:
					process[i] = 2;
					break;
				default:
					process[i] = 1;
					break;
				}
			}
			break;
		}
		case 1:
		{
			diameter[i] = diameterU;
			sDev[i] = sdevU;
			iDev[i] = idevU;
			cs[i] = csU;
			if (vi.IsRGB())
			{
				process[i] = 3;
				width[i] = vi.width;
				height[i] = vi.height;
				if (diameter[i] > width[i] || diameter[i] > height[i])
					env->ThrowError("vsTBilateral: diameter must be less than the dimensions of the plane.");
			}
			else
			{
				width[i] = vi.width >> vi.GetPlaneWidthSubsampling(PLANAR_U);
				height[i] = vi.height >> vi.GetPlaneHeightSubsampling(PLANAR_U);
				if (diameter[i] > width[i] || diameter[i] > height[i])
					env->ThrowError("vsTBilateral: diameterU must be less than the dimensions of plane U.");

				switch (u_)
				{
				case 3:
					process[i] = 3;
					break;
				case 2:
					process[i] = 2;
					break;
				default:
					process[i] = 1;
					break;
				}
			}
			break;
		}
		default:
		{
			diameter[i] = diameterV;
			sDev[i] = sdevV;
			iDev[i] = idevV;
			cs[i] = csV;
			if (vi.IsRGB())
			{
				process[i] = 3;
				width[i] = vi.width;
				height[i] = vi.height;
				if (diameter[i] > width[i] || diameter[i] > height[i])
					env->ThrowError("vsTBilateral: diameter must be less than the dimensions of the plane.");
			}
			else
			{
				width[i] = vi.width >> vi.GetPlaneWidthSubsampling(PLANAR_V);
				height[i] = vi.height >> vi.GetPlaneHeightSubsampling(PLANAR_V);
				if (diameter[i] > width[i] || diameter[i] > height[i])
					env->ThrowError("vsTBilateral: diameterV must be less than the dimensions of plane V.");

				switch (v_)
				{
				case 3:
					process[i] = 3;
					break;
				case 2:
					process[i] = 2;
					break;
				default:
					process[i] = 1;
					break;
				}
			}
			break;
		}
		}

		if (diameter[i] <= 1 || diameter[i] % 2 == 0)
			env->ThrowError("vsTBilateral: diameterY/U/V must be an odd number greater than 1.");
		if (sDev[i] < 0)
			env->ThrowError("vsTBilateral: sdevY/U/V must be at least 0.");
		if (iDev[i] < 0)
			env->ThrowError("vsTBilateral: idevY/U/V must be at least 0.");
		if (cs[i] < 0)
			env->ThrowError("vsTBilateral: csY/U/V must be at least 0.");

		int window = diameter[i] * diameter[i];
		int radius = diameter[i] >> 1;

		spatialWeights[i] = reinterpret_cast<double*>(_aligned_malloc(window * sizeof(double), 16));
		double* disTable = reinterpret_cast<double*>(_aligned_malloc(window * sizeof(double), 16));

		for (int b = 0, y = -radius; y <= radius; ++y)
		{
			int temp = y * y;
			for (int x = -radius; x <= radius; ++x)
				disTable[b++] = sqrt((double)(temp + static_cast<int64_t>(x) * x));
		}

		for (int x = 0; x < window; ++x)
			spatialWeights[i][x] = kernelValue(disTable[x], sDev[i], kerns_);
		spatialWeights[i][radius * diameter[i] + radius] *= cs[i];

		int diff_size = pixel_max * (d2_ ? 4 : 2) + 1;

		diffWeights[i] = reinterpret_cast<double*>(_aligned_malloc(diff_size * sizeof(double), 16));

		iDev[i] = iDev[i] * pixel_max / 255;

		for (int x = 0; x <= diff_size / 2; ++x)
			diffWeights[i][diff_size / 2 + x] = diffWeights[i][diff_size / 2 - x] = kernelValue(x / (d2_ ? 2.0 : 1.0), iDev[i], kerni_);

		_aligned_free(disTable);
	}

	has_at_least_v8 = true;
	try { env->CheckVersion(8); }
	catch (const AvisynthError&) { has_at_least_v8 = false; }
}

vsTBilateral::~vsTBilateral()
{
	const int planecount = std::min(vi.NumComponents(), 3);
	for (int i = 0; i < planecount; i++)
	{
		_aligned_free(spatialWeights[i]);
		_aligned_free(diffWeights[i]);
	}
}

PVideoFrame __stdcall vsTBilateral::GetFrame(int n, IScriptEnvironment* env)
{
	PVideoFrame src = child->GetFrame(n, env);
	PVideoFrame pp = ppclip_ ? ppclip_->GetFrame(n, env) : src;
	PVideoFrame dst = has_at_least_v8 ? env->NewVideoFrameP(vi, &src) : env->NewVideoFrame(vi);

	int planes_y[4] = { PLANAR_Y, PLANAR_U, PLANAR_V, PLANAR_A };
	int planes_r[4] = { PLANAR_G, PLANAR_B, PLANAR_R, PLANAR_A };
	const int* current_planes = (vi.IsYUV() || vi.IsYUVA()) ? planes_y : planes_r;
	const int planecount = std::min(vi.NumComponents(), 3);
	for (int i = 0; i < planecount; i++)
	{
		const int plane = current_planes[i];

		switch (process[i])
		{
		case 3:
		{
			int bits = vi.BitsPerComponent();
			if (d2_)
			{
				switch (resType_)
				{
				case 0:
				{
					switch (bits)
					{
					case 8:
						ProcessFrameD2_Mean<uint8_t>(src->GetReadPtr(plane), pp->GetReadPtr(plane), dst->GetWritePtr(plane), src->GetPitch(plane), pp->GetPitch(plane), dst->GetPitch(plane), src->GetRowSize(plane) / vi.ComponentSize(), src->GetHeight(plane), i, env);
						break;
					default:
						ProcessFrameD2_Mean<uint16_t>(src->GetReadPtr(plane), pp->GetReadPtr(plane), dst->GetWritePtr(plane), src->GetPitch(plane), pp->GetPitch(plane), dst->GetPitch(plane), src->GetRowSize(plane) / vi.ComponentSize(), src->GetHeight(plane), i, env);
						break;
					}
					break;
				}
				case 3:
				{
					switch (bits)
					{
					case 8:
						ProcessFrameD2_MLR<uint8_t>(src->GetReadPtr(plane), pp->GetReadPtr(plane), dst->GetWritePtr(plane), src->GetPitch(plane), pp->GetPitch(plane), dst->GetPitch(plane), src->GetRowSize(plane) / vi.ComponentSize(), src->GetHeight(plane), i, env);
						break;
					default:
						ProcessFrameD2_MLR<uint16_t>(src->GetReadPtr(plane), pp->GetReadPtr(plane), dst->GetWritePtr(plane), src->GetPitch(plane), pp->GetPitch(plane), dst->GetPitch(plane), src->GetRowSize(plane) / vi.ComponentSize(), src->GetHeight(plane), i, env);
						break;
					}
					break;
				}
				default:
				{
					switch (bits)
					{
					case 8:
						ProcessFrameD2_Med<uint8_t>(src->GetReadPtr(plane), pp->GetReadPtr(plane), dst->GetWritePtr(plane), src->GetPitch(plane), pp->GetPitch(plane), dst->GetPitch(plane), src->GetRowSize(plane) / vi.ComponentSize(), src->GetHeight(plane), i, env);
						break;
					default:
						ProcessFrameD2_Med<uint16_t>(src->GetReadPtr(plane), pp->GetReadPtr(plane), dst->GetWritePtr(plane), src->GetPitch(plane), pp->GetPitch(plane), dst->GetPitch(plane), src->GetRowSize(plane) / vi.ComponentSize(), src->GetHeight(plane), i, env);
						break;
					}
					break;
				}
				}
			}
			else
			{
				switch (resType_)
				{
				case 0:
				{
					switch (bits)
					{
					case 8:
						ProcessFrameD1_Mean<uint8_t>(src->GetReadPtr(plane), pp->GetReadPtr(plane), dst->GetWritePtr(plane), src->GetPitch(plane), pp->GetPitch(plane), dst->GetPitch(plane), src->GetRowSize(plane) / vi.ComponentSize(), src->GetHeight(plane), i, env);
						break;
					default:
						ProcessFrameD1_Mean<uint16_t>(src->GetReadPtr(plane), pp->GetReadPtr(plane), dst->GetWritePtr(plane), src->GetPitch(plane), pp->GetPitch(plane), dst->GetPitch(plane), src->GetRowSize(plane) / vi.ComponentSize(), src->GetHeight(plane), i, env);
						break;
					}
					break;
				}
				case 3:
				{
					switch (bits)
					{
					case 8:
						ProcessFrameD1_MLR<uint8_t>(src->GetReadPtr(plane), pp->GetReadPtr(plane), dst->GetWritePtr(plane), src->GetPitch(plane), pp->GetPitch(plane), dst->GetPitch(plane), src->GetRowSize(plane) / vi.ComponentSize(), src->GetHeight(plane), i, env);
						break;
					default:
						ProcessFrameD1_MLR<uint16_t>(src->GetReadPtr(plane), pp->GetReadPtr(plane), dst->GetWritePtr(plane), src->GetPitch(plane), pp->GetPitch(plane), dst->GetPitch(plane), src->GetRowSize(plane) / vi.ComponentSize(), src->GetHeight(plane), i, env);
						break;
					}
					break;
				}
				default:
				{
					switch (bits)
					{
					case 8:
						ProcessFrameD1_Med<uint8_t>(src->GetReadPtr(plane), pp->GetReadPtr(plane), dst->GetWritePtr(plane), src->GetPitch(plane), pp->GetPitch(plane), dst->GetPitch(plane), src->GetRowSize(plane) / vi.ComponentSize(), src->GetHeight(plane), i, env);
						break;
					default:
						ProcessFrameD1_Med<uint16_t>(src->GetReadPtr(plane), pp->GetReadPtr(plane), dst->GetWritePtr(plane), src->GetPitch(plane), pp->GetPitch(plane), dst->GetPitch(plane), src->GetRowSize(plane) / vi.ComponentSize(), src->GetHeight(plane), i, env);
						break;
					}
					break;
				}
				}
			}
			break;
		}
		case 2:
			env->BitBlt(dst->GetWritePtr(plane), dst->GetPitch(plane), src->GetReadPtr(plane), src->GetPitch(plane), src->GetRowSize(plane), src->GetHeight(plane));
			break;
		}
	}

	return dst;
}

AVSValue __cdecl Create_vsTBilateral(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	return new vsTBilateral(
		args[0].AsClip(),
		args[1].AsClip(),
		args[2].AsInt(5),
		args[3].AsInt(5),
		args[4].AsInt(5),
		args[5].AsFloat(1.4f),
		args[6].AsFloat(1.4f),
		args[7].AsFloat(1.4f),
		args[8].AsFloat(7.0f),
		args[9].AsFloat(7.0f),
		args[10].AsFloat(7.0f),
		args[11].AsFloat(1.0f),
		args[12].AsFloat(1.0f),
		args[13].AsFloat(1.0f),
		args[14].AsBool(false),
		args[15].AsInt(2),
		args[16].AsInt(2),
		args[17].AsInt(0),
		args[18].AsInt(3),
		args[19].AsInt(3),
		args[20].AsInt(3),
		env);
}

const AVS_Linkage* AVS_linkage;

extern "C" __declspec(dllexport)
const char* __stdcall AvisynthPluginInit3(IScriptEnvironment * env, const AVS_Linkage* const vectors)
{
	AVS_linkage = vectors;

	env->AddFunction("vsTBilateral", "c[ppclip]c[diameterY]i[diameterU]i[diameterV]i[sdevY]f[sdevU]f[sdevV]f[idevY]f[idevU]f[idevV]f[csY]f[csU]f[csV]f[d2]b[kerns]i[kerni]i[restype]i[y]i[u]i[v]i", Create_vsTBilateral, 0);

	return "vsTBilateral";
}

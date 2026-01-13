/*
 * MeshExport - Exports obj and x3d files with color from processing with beginRecord and endRecord
 * by Jesse Louis-Rosenberg 
 * http://n-e-r-v-o-u-s.com
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * <p/>
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * <p/>
 * You should have received a copy of the GNU Lesser General
 * Public License; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place,
 * Suite 330, Boston, MA  02111-1307  USA
 */

package nervoussystem.obj;

import java.io.*;
import java.util.HashMap;
import java.awt.*;
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;
import processing.core.*;

public class MeshExport extends PGraphics {
  File file;
  String filenameSimple;
  PrintWriter writer;
  float[][] pts;
  int[][] lines;
  int[][] faces;
  int[][] faceColors;
  int[] colors;
  int lineCount;
  int faceCount;
  int objectMode = 0;
  HashMap<String,Integer> ptMap;
  
  static protected final int MATRIX_STACK_DEPTH = 32;
  int DEFAULT_VERTICES = 4096;
  int VERTEX_FIELD_COUNT = 3;
  int vertices[];
  
  int numTriangles = 0;
  int numQuads = 0;
  int ptCount = 0;
  static public int TRIANGLE_RES = 7;
  static public int RECT_RES = TRIANGLE_RES+5;
  
  boolean drawBegan = false;
  //make transform function work
  protected int stackDepth;
  protected float[][] matrixStack = new float[MATRIX_STACK_DEPTH][16];
  PMatrix3D transformMatrix;
  
  //color
  boolean colorFlag = false;
  private static BufferedImage textureImg = null;
  private static Graphics2D g2d = null;
  private static int textureWidth = 0;
  private static int textureHeight = 0;
  
  public MeshExport() {
    vertices = new int[DEFAULT_VERTICES];
	stackDepth = 0;
	transformMatrix = new PMatrix3D();
  }

  public void setPath(String path) {
    this.path = path;
    if (path != null) {
      file = new File(path);
      if (!file.isAbsolute()) file = null;
    }
    if (file == null) {
      throw new RuntimeException("OBJExport requires an absolute path " +
        "for the location of the output file.");
    }
	filenameSimple = file.getName();
	int dotPos = filenameSimple.lastIndexOf(".");
	if(dotPos > -1) {
		filenameSimple = filenameSimple.substring(0,dotPos);
	}
  }

  protected void allocate() {
  }

  public void dispose() {
    writer.flush();
    writer.close();
    writer = null;
	ptMap.clear();
	if(textureImg != null) {
		textureImg = null;
		g2d = null;
	}
  }

  public boolean displayable() {
    return false;  // just in case someone wants to use this on its own
  }

  public void beginDraw() {
    // have to create file object here, because the name isn't yet
    // available in allocate()
    // Processing 4.x compatibility: removed defaultSettings() call
    // as it causes NullPointerException when parent is not yet set.
    // Initialize only what we need for this renderer.
    colorMode(RGB, 255);
    fillColor = 0xFFFFFFFF;  // white fill by default
    fill = true;
    shape = 0;

    if (writer == null) {
      try {
        writer = new PrintWriter(new FileWriter(file));
      } 
      catch (IOException e) {
        throw new RuntimeException(e);
      }
      pts = new float[4096][VERTEX_FIELD_COUNT];
      lines = new int[4096][];
      faces = new int[4096][];
      ptMap = new HashMap<String,Integer>();
	  if(colorFlag) {
		colors = new int[4096];
		faceColors = new int[4096][];
	  }
    }
    lineCount = 0;
    faceCount = 0;
    vertexCount = 0;
	ptCount = 0;
    numTriangles = 0;
    numQuads = 0;
  }

  public void endDraw() {
    //write vertices and initialize ptMap
	writeHeader();
	if(colorFlag) {
		colorExport();
	}
    writeVertices();
    writeLines();
    writeFaces();
	writeFooter();
	drawBegan = false;
  }
    
  private void colorExport() {
	int numRects = PApplet.ceil(numTriangles/2.0f) + numQuads;
    int textureSize = PApplet.ceil(PApplet.sqrt(numRects))*RECT_RES;

	if(textureSize > 1024) {
		showWarning("Generating texture... this might take a while");
		textureSize = PApplet.ceil(textureSize/1024.0f)*1024;
	}

	// Use BufferedImage directly to avoid pixelDensity issues
	textureWidth = textureSize;
	textureHeight = textureSize;
	textureImg = new BufferedImage(textureSize, textureSize, BufferedImage.TYPE_INT_RGB);
	g2d = textureImg.createGraphics();
	g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF);
	g2d.setColor(Color.BLACK);
	g2d.fillRect(0, 0, textureSize, textureSize);

	generateTexture();
  }
      
  private void generateTexture() {
	int c;
	int currX = 0, currY = 0;
	int[] f;
	int[] fc;
	boolean upper = true;

	for(int i=0;i<faceCount;++i) {
		f = faces[i];
		fc = faceColors[i];
		if(f.length > 4) showWarning("Faces with more than 4 sides cannot be exported with color");
		else if(f.length == 3) {
			//draw triangle
			if(upper) {
				int x0 = currX+1, y0 = currY+1;
				int x1 = currX+1+TRIANGLE_RES, y1 = currY+1;
				int x2 = currX+1, y2 = currY+TRIANGLE_RES+1;
				int[] xPoints = {x0, x1, x2};
				int[] yPoints = {y0, y1, y2};
				// Fill with gradient approximation - use vertex colors
				drawColoredTriangle(x0, y0, fc[0], x1, y1, fc[1], x2, y2, fc[2]);
				writeTexCoord((x0+.5f)/textureWidth, (1.0f-(y0+.5f)/textureHeight));
				writeTexCoord((x1+.5f)/textureWidth, (1.0f-(y1+.5f)/textureHeight));
				writeTexCoord((x2+.5f)/textureWidth, (1.0f-(y2+.5f)/textureHeight));
			} else {
				int x0 = currX+3+TRIANGLE_RES, y0 = currY+3;
				int x1 = currX+TRIANGLE_RES+3, y1 = currY+TRIANGLE_RES+3;
				int x2 = currX+3, y2 = currY+TRIANGLE_RES+3;
				drawColoredTriangle(x0, y0, fc[0], x1, y1, fc[1], x2, y2, fc[2]);
				writeTexCoord((x0+.5f)/textureWidth, (1.0f-(y0+.5f)/textureHeight));
				writeTexCoord((x1+.5f)/textureWidth, (1.0f-(y1+.5f)/textureHeight));
				writeTexCoord((x2+.5f)/textureWidth, (1.0f-(y2+.5f)/textureHeight));
				currX += RECT_RES;
			}
			upper = !upper;
		} else if(f.length == 4) {
			int x0 = currX+1, y0 = currY+1;
			int x1 = currX+1+TRIANGLE_RES, y1 = currY+1;
			int x2 = currX+TRIANGLE_RES+1, y2 = currY+TRIANGLE_RES+1;
			int x3 = currX+1, y3 = currY+TRIANGLE_RES+1;
			drawColoredQuad(x0, y0, fc[0], x1, y1, fc[1], x2, y2, fc[2], x3, y3, fc[3]);
			writeTexCoord((x0+.5f)/textureWidth, (1.0f-(y0+.5f)/textureHeight));
			writeTexCoord((x1+.5f)/textureWidth, (1.0f-(y1+.5f)/textureHeight));
			writeTexCoord((x2+.5f)/textureWidth, (1.0f-(y2+.5f)/textureHeight));
			writeTexCoord((x3+.5f)/textureWidth, (1.0f-(y3+.5f)/textureHeight));
			currX += RECT_RES;
			upper = true;
		}
		if(currX+RECT_RES > textureWidth) {
			currX = 0;
			currY += RECT_RES;
		}
	}

	g2d.dispose();

	// Save texture
	try {
		ImageIO.write(textureImg, "PNG", new File(file.getParent() + File.separator + filenameSimple + ".png"));
	} catch (IOException e) {
		showWarning("Failed to save texture: " + e.getMessage());
	}
  }

  // Draw a triangle with vertex color interpolation using barycentric coordinates
  private void drawColoredTriangle(int x0, int y0, int c0, int x1, int y1, int c1, int x2, int y2, int c2) {
	// Calculate bounding box
	int minX = Math.min(x0, Math.min(x1, x2));
	int maxX = Math.max(x0, Math.max(x1, x2));
	int minY = Math.min(y0, Math.min(y1, y2));
	int maxY = Math.max(y0, Math.max(y1, y2));

	// Precompute denominator for barycentric coordinates
	float denom = (float)((y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2));
	if (Math.abs(denom) < 0.0001f) {
		// Degenerate triangle, use average color
		int r = (((c0 >> 16) & 0xFF) + ((c1 >> 16) & 0xFF) + ((c2 >> 16) & 0xFF)) / 3;
		int g = (((c0 >> 8) & 0xFF) + ((c1 >> 8) & 0xFF) + ((c2 >> 8) & 0xFF)) / 3;
		int b = ((c0 & 0xFF) + (c1 & 0xFF) + (c2 & 0xFF)) / 3;
		g2d.setColor(new Color(r, g, b));
		g2d.fillPolygon(new int[]{x0, x1, x2}, new int[]{y0, y1, y2}, 3);
		return;
	}

	// Extract RGB components from vertex colors
	int r0 = (c0 >> 16) & 0xFF, g0 = (c0 >> 8) & 0xFF, b0 = c0 & 0xFF;
	int r1 = (c1 >> 16) & 0xFF, g1 = (c1 >> 8) & 0xFF, b1 = c1 & 0xFF;
	int r2 = (c2 >> 16) & 0xFF, g2 = (c2 >> 8) & 0xFF, b2 = c2 & 0xFF;

	// Iterate over bounding box and interpolate colors
	for (int py = minY; py <= maxY; py++) {
		for (int px = minX; px <= maxX; px++) {
			// Bounds check
			if (px < 0 || px >= textureWidth || py < 0 || py >= textureHeight) {
				continue;
			}

			// Calculate barycentric coordinates
			float lambda0 = ((y1 - y2) * (px - x2) + (x2 - x1) * (py - y2)) / denom;
			float lambda1 = ((y2 - y0) * (px - x2) + (x0 - x2) * (py - y2)) / denom;
			float lambda2 = 1.0f - lambda0 - lambda1;

			// Check if point is inside triangle (with small epsilon for edge cases)
			if (lambda0 >= -0.001f && lambda1 >= -0.001f && lambda2 >= -0.001f) {
				// Clamp lambdas to valid range
				lambda0 = Math.max(0, Math.min(1, lambda0));
				lambda1 = Math.max(0, Math.min(1, lambda1));
				lambda2 = Math.max(0, Math.min(1, lambda2));

				// Normalize
				float sum = lambda0 + lambda1 + lambda2;
				lambda0 /= sum;
				lambda1 /= sum;
				lambda2 /= sum;

				// Interpolate color
				int r = (int)(lambda0 * r0 + lambda1 * r1 + lambda2 * r2);
				int g = (int)(lambda0 * g0 + lambda1 * g1 + lambda2 * g2);
				int b = (int)(lambda0 * b0 + lambda1 * b1 + lambda2 * b2);

				// Clamp to valid range
				r = Math.max(0, Math.min(255, r));
				g = Math.max(0, Math.min(255, g));
				b = Math.max(0, Math.min(255, b));

				textureImg.setRGB(px, py, (r << 16) | (g << 8) | b);
			}
		}
	}

	// Draw border to fill edge pixels
	drawTriangleBorder(x0, y0, c0, x1, y1, c1, x2, y2, c2);
  }

  // Draw triangle border with interpolated colors along edges
  private void drawTriangleBorder(int x0, int y0, int c0, int x1, int y1, int c1, int x2, int y2, int c2) {
	drawInterpolatedLine(x0, y0, c0, x1, y1, c1);
	drawInterpolatedLine(x1, y1, c1, x2, y2, c2);
	drawInterpolatedLine(x2, y2, c2, x0, y0, c0);
  }

  // Draw a line with color interpolation between two endpoints
  private void drawInterpolatedLine(int x0, int y0, int c0, int x1, int y1, int c1) {
	int dx = Math.abs(x1 - x0);
	int dy = Math.abs(y1 - y0);
	int steps = Math.max(dx, dy);
	if (steps == 0) {
		if (x0 >= 0 && x0 < textureWidth && y0 >= 0 && y0 < textureHeight) {
			textureImg.setRGB(x0, y0, c0 & 0xFFFFFF);
		}
		return;
	}

	int r0 = (c0 >> 16) & 0xFF, g0 = (c0 >> 8) & 0xFF, b0 = c0 & 0xFF;
	int r1 = (c1 >> 16) & 0xFF, g1 = (c1 >> 8) & 0xFF, b1 = c1 & 0xFF;

	for (int i = 0; i <= steps; i++) {
		float t = (float) i / steps;
		int px = (int)(x0 + t * (x1 - x0));
		int py = (int)(y0 + t * (y1 - y0));

		int r = (int)(r0 + t * (r1 - r0));
		int g = (int)(g0 + t * (g1 - g0));
		int b = (int)(b0 + t * (b1 - b0));

		if (px >= 0 && px < textureWidth && py >= 0 && py < textureHeight) {
			textureImg.setRGB(px, py, (r << 16) | (g << 8) | b);
		}
	}
  }

  // Draw a quad with vertex color interpolation (using bilinear interpolation)
  private void drawColoredQuad(int x0, int y0, int c0, int x1, int y1, int c1, int x2, int y2, int c2, int x3, int y3, int c3) {
	// Calculate bounding box
	int minX = Math.min(Math.min(x0, x1), Math.min(x2, x3));
	int maxX = Math.max(Math.max(x0, x1), Math.max(x2, x3));
	int minY = Math.min(Math.min(y0, y1), Math.min(y2, y3));
	int maxY = Math.max(Math.max(y0, y1), Math.max(y2, y3));

	// Extract RGB components from vertex colors
	int r0 = (c0 >> 16) & 0xFF, g0 = (c0 >> 8) & 0xFF, b0 = c0 & 0xFF;
	int r1 = (c1 >> 16) & 0xFF, g1 = (c1 >> 8) & 0xFF, b1 = c1 & 0xFF;
	int r2 = (c2 >> 16) & 0xFF, g2 = (c2 >> 8) & 0xFF, b2 = c2 & 0xFF;
	int r3 = (c3 >> 16) & 0xFF, g3 = (c3 >> 8) & 0xFF, b3 = c3 & 0xFF;

	// For each pixel in bounding box, check if inside quad and interpolate color
	for (int py = minY; py <= maxY; py++) {
		for (int px = minX; px <= maxX; px++) {
			// Bounds check
			if (px < 0 || px >= textureWidth || py < 0 || py >= textureHeight) {
				continue;
			}

			// Check if point is inside the quad using cross product test
			if (!isInsideQuad(px, py, x0, y0, x1, y1, x2, y2, x3, y3)) {
				continue;
			}

			// Use bilinear interpolation based on normalized position
			// Quad vertices: 0(top-left), 1(top-right), 2(bottom-right), 3(bottom-left)
			float u = (float)(px - minX) / Math.max(1, maxX - minX);
			float v = (float)(py - minY) / Math.max(1, maxY - minY);

			// Bilinear interpolation
			// top edge: lerp between c0 and c1
			// bottom edge: lerp between c3 and c2
			// then lerp between top and bottom
			float rTop = r0 + u * (r1 - r0);
			float gTop = g0 + u * (g1 - g0);
			float bTop = b0 + u * (b1 - b0);

			float rBottom = r3 + u * (r2 - r3);
			float gBottom = g3 + u * (g2 - g3);
			float bBottom = b3 + u * (b2 - b3);

			int r = (int)(rTop + v * (rBottom - rTop));
			int g = (int)(gTop + v * (gBottom - gTop));
			int b = (int)(bTop + v * (bBottom - bTop));

			// Clamp to valid range
			r = Math.max(0, Math.min(255, r));
			g = Math.max(0, Math.min(255, g));
			b = Math.max(0, Math.min(255, b));

			textureImg.setRGB(px, py, (r << 16) | (g << 8) | b);
		}
	}

	// Draw border to fill edge pixels
	drawQuadBorder(x0, y0, c0, x1, y1, c1, x2, y2, c2, x3, y3, c3);
  }

  // Check if a point is inside a quadrilateral using cross product test
  private boolean isInsideQuad(int px, int py, int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3) {
	// Check if point is on the same side of all edges
	float d0 = crossProduct2D(x1 - x0, y1 - y0, px - x0, py - y0);
	float d1 = crossProduct2D(x2 - x1, y2 - y1, px - x1, py - y1);
	float d2 = crossProduct2D(x3 - x2, y3 - y2, px - x2, py - y2);
	float d3 = crossProduct2D(x0 - x3, y0 - y3, px - x3, py - y3);

	boolean hasNeg = (d0 < 0) || (d1 < 0) || (d2 < 0) || (d3 < 0);
	boolean hasPos = (d0 > 0) || (d1 > 0) || (d2 > 0) || (d3 > 0);

	// If all same sign (or zero), point is inside
	return !(hasNeg && hasPos);
  }

  // 2D cross product
  private float crossProduct2D(float ax, float ay, float bx, float by) {
	return ax * by - ay * bx;
  }

  // Draw quad border with interpolated colors along edges
  private void drawQuadBorder(int x0, int y0, int c0, int x1, int y1, int c1, int x2, int y2, int c2, int x3, int y3, int c3) {
	drawInterpolatedLine(x0, y0, c0, x1, y1, c1);
	drawInterpolatedLine(x1, y1, c1, x2, y2, c2);
	drawInterpolatedLine(x2, y2, c2, x3, y3, c3);
	drawInterpolatedLine(x3, y3, c3, x0, y0, c0);
  }
  
  public void beginShape(int kind) {
    shape = kind;
    vertexCount = 0;
  }

  public void vertex(float x, float y) {
    vertex(x,y,0);
  }
  
  private static float tempVertex[] = new float[3];
  
  float[] vertex = new float[3];
  public void vertex(float x, float y,float z) {	
    if(vertexCount >= vertices.length) {
		int newVertices[] = new int[vertices.length*2];
		System.arraycopy(vertices,0,newVertices,0,vertices.length);
		vertices = newVertices;
    }
   // float vertex[] = vertices[vertexCount];
    vertex[X] = x;
    vertex[Y] = y;
    vertex[Z] = z;
	transformMatrix.mult(vertex, tempVertex);
	x = tempVertex[0];
	y = tempVertex[1];
	z = tempVertex[2];
	//does not account for floating point error or tolerance
    if(!ptMap.containsKey(x+"_"+y+"_"+z)) {
      if(ptCount >= pts.length) {
		float newPts[][] = new float[pts.length*2][];
		System.arraycopy(pts,0,newPts,0,pts.length);
		pts = newPts;
      }
	  //might need to separating position and color so faces can have different colors
	  //the plan: make every call of fill add a new color, have a separate uv index for the faces
	  pts[ptCount] = new float[] {x,y,z};
	  ptCount++;
	  ptMap.put(x+"_"+y+"_"+z,new Integer(ptMap.size()+1));
	  vertices[vertexCount] = ptCount;
    } else {
		vertices[vertexCount] = ptMap.get(x+"_"+y+"_"+z).intValue();
	}
	//color
	if(colorFlag) {
		if(vertexCount >= colors.length) {
			int newColors[] = new int[colors.length*2];
			System.arraycopy(colors,0,newColors,0,colors.length);
			colors = newColors;
		}
		colors[vertexCount] = fillColor;
	}

    vertexCount++;
  }

  public void endShape(int mode) {
    //if(stroke) endShapeStroke(mode);
    //if(fill) endShapeFill(mode);
    endShapeFill(mode);
 }

  
  public void endShapeFill(int mode) {
      switch(shape) {
      case TRIANGLES:
        {
        int stop = vertexCount-2;
          for (int i = 0; i < stop; i += 3) {
            int[] f = new int[3];
            f[0] = vertices[i];
            f[1] = vertices[i+1];
            f[2] = vertices[i+2];
            addFace(f);
			if(colorFlag) {
				int[] fc = new int[3];
				fc[0] = colors[i];
				fc[1] = colors[i+1];
				fc[2] = colors[i+2];
				addFaceColor(fc,faceCount-1);
			}
          }
        }
        break;
      case TRIANGLE_STRIP:
      {
          int stop = vertexCount - 2;
          for (int i = 0; i < stop; i++) {
            // have to switch between clockwise/counter-clockwise
            // otherwise the feller is backwards and renderer won't draw
            if ((i % 2) == 0) {
              int[] f = new int[3];
              f[0] = vertices[i];
              f[1] = vertices[i+2];
              f[2] = vertices[i+1];
              addFace(f);
			  if(colorFlag) {
				int[] fc = new int[3];
				fc[0] = colors[i];
				fc[1] = colors[i+2];
				fc[2] = colors[i+1];
				addFaceColor(fc,faceCount-1);
			  }
            } else {
              int[] f = new int[3];
              f[0] = vertices[i];
              f[1] = vertices[i+1];
              f[2] = vertices[i+2];
              addFace(f);
			  if(colorFlag) {
				int[] fc = new int[3];
				fc[0] = colors[i];
				fc[1] = colors[i+1];
				fc[2] = colors[i+2];
				addFaceColor(fc, faceCount-1);
			  }
            }
          }
      }
      break;
      case POLYGON:
      {
        int[] f;
        boolean closed = vertices[0]!=vertices[vertexCount-1];
        if(closed) {
         f = new int[vertexCount];
        } else {
         f = new int[vertexCount-1];
        }
        int end = vertexCount;
        if(!closed) end--;
        for(int i=0;i<end;++i) {
          f[i] = vertices[i];
        }
        addFace(f);
		if(colorFlag) {
			if(end <= 4) {
				int[] fc = new int[end];
				for(int i=0;i<end;++i) {
				  fc[i] = colors[i];
				}
				addFaceColor(fc,faceCount-1);
			} else {
				//dummy so faces and faceColors match (stupid)
				addFaceColor(new int[0],faceCount-1);
			}
		}
      }
      break;
      case QUADS:
      {
        int stop = vertexCount-3;
        for (int i = 0; i < stop; i += 4) {
            int[] f = new int[4];
            f[0] = vertices[i];
            f[1] = vertices[i+1];
            f[2] = vertices[i+2];
            f[3] = vertices[i+3];
            addFace(f);
			if(colorFlag) {
				int[] fc = new int[4];
				fc[0] = colors[i];
				fc[1] = colors[i+1];
				fc[2] = colors[i+2];
				fc[3] = colors[i+3];
				addFaceColor(fc,faceCount-1);
			}
        }
      }
      break;

      case QUAD_STRIP:
      {
        int stop = vertexCount-3;
        for (int i = 0; i < stop; i += 2) {
            int[] f = new int[4];
            f[0] = vertices[i];
            f[1] = vertices[i+1];
            f[3] = vertices[i+2];
            f[2] = vertices[i+3];
            addFace(f);        
			if(colorFlag) {
				int[] fc = new int[4];
				fc[0] = colors[i];
				fc[1] = colors[i+1];
				fc[2] = colors[i+2];
				fc[3] = colors[i+3];
				addFaceColor(fc,faceCount-1);
			}
		}
      }
      break;
      case TRIANGLE_FAN:
      {
        int stop = vertexCount - 1;
        for (int i = 1; i < stop; i++) {
			int f[] = new int[3];
			f[0] = vertices[0];
			f[1] = vertices[i];
			f[2] = vertices[i+1];
			addFace(f);
			if(colorFlag) {
				int[] fc = new int[3];
				fc[0] = colors[0];
				fc[1] = colors[i];
				fc[2] = colors[i+1];
				addFaceColor(fc,faceCount-1);
			}
		}
      }
      break;
    }
  }
  
  //unused as of this version
  public void endShapeStroke(int mode) {
  }
  
  private void addFace(int[] f) {
   if(faceCount >= faces.length) {
    int newfaces[][] = new int[faces.length*2][];
    System.arraycopy(faces,0,newfaces,0,faces.length);
    faces = newfaces;
   }
   if(f.length == 3) numTriangles++;
   else if(f.length == 4) numQuads++;
   faces[faceCount++] = f;
  }
  
  private void addFaceColor(int[] f, int pos) {
   if(pos >= faceColors.length) {
    int newfaces[][] = new int[faceColors.length*2][];
    System.arraycopy(faceColors,0,newfaces,0,faceColors.length);
    faceColors = newfaces;
   }
   faceColors[pos] = f;
  }
  
  private void addLine(int[] l) {
   if(lineCount >= lines.length) {
    int newLines[][] = new int[lines.length*2][];
    System.arraycopy(lines,0,newLines,0,lines.length);
    lines = newLines;
   }
   lines[lineCount++] = l;
  }  
  
  //taken from PGraphicsOpenGL
  @Override
  public void pushMatrix() {
    if (stackDepth == MATRIX_STACK_DEPTH) {
      throw new RuntimeException(ERROR_PUSHMATRIX_OVERFLOW);
    }
    transformMatrix.get(matrixStack[stackDepth]);
    stackDepth++;
  }

  @Override
  public void popMatrix() {
    if (stackDepth == 0) {
      throw new RuntimeException(ERROR_PUSHMATRIX_UNDERFLOW);
    }
    stackDepth--;
    transformMatrix.set(matrixStack[stackDepth]);
  }
  
  @Override
  public void translate(float tx, float ty) {
    translateImpl(tx, ty, 0);
  }

  @Override
  public void translate(float tx, float ty, float tz) {
    translateImpl(tx, ty, tz);
  }

  protected void translateImpl(float tx, float ty, float tz) {
    transformMatrix.translate(tx, ty, tz);
  }
  
  /**
   * Two dimensional rotation. Same as rotateZ (this is identical to a 3D
   * rotation along the z-axis) but included for clarity -- it'd be weird for
   * people drawing 2D graphics to be using rotateZ. And they might kick our a--
   * for the confusion.
   */
  @Override
  public void rotate(float angle) {
    rotateImpl(angle, 0, 0, 1);
  }


  @Override
  public void rotateX(float angle) {
    rotateImpl(angle, 1, 0, 0);
  }


  @Override
  public void rotateY(float angle) {
    rotateImpl(angle, 0, 1, 0);
  }


  @Override
  public void rotateZ(float angle) {
    rotateImpl(angle, 0, 0, 1);
  }


  /**
   * Rotate around an arbitrary vector, similar to glRotate(), except that it
   * takes radians (instead of degrees).
   */
  @Override
  public void rotate(float angle, float v0, float v1, float v2) {
    rotateImpl(angle, v0, v1, v2);
  }


  protected void rotateImpl(float angle, float v0, float v1, float v2) {
    float norm2 = v0 * v0 + v1 * v1 + v2 * v2;
    if (zero(norm2)) {
      // The vector is zero, cannot apply rotation.
      return;
    }

    if (diff(norm2, 1)) {
      // The rotation vector is not normalized.
      float norm = PApplet.sqrt(norm2);
      v0 /= norm;
      v1 /= norm;
      v2 /= norm;
    }

    transformMatrix.rotate(angle, v0, v1, v2);
  }
  
  /**
   * Same as scale(s, s, s).
   */
  @Override
  public void scale(float s) {
    scaleImpl(s, s, s);
  }


  /**
   * Same as scale(sx, sy, 1).
   */
  @Override
  public void scale(float sx, float sy) {
    scaleImpl(sx, sy, 1);
  }


  /**
   * Scale in three dimensions.
   */
  @Override
  public void scale(float sx, float sy, float sz) {
    scaleImpl(sx, sy, sz);
  }

  /**
   * Scale in three dimensions.
   */
  protected void scaleImpl(float sx, float sy, float sz) {
    transformMatrix.scale(sx, sy, sz);
  }
  
  @Override
  public void shearX(float angle) {
    float t = (float) Math.tan(angle);
    applyMatrixImpl(1, t, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1);
  }


  @Override
  public void shearY(float angle) {
    float t = (float) Math.tan(angle);
    applyMatrixImpl(1, 0, 0, 0,
                    t, 1, 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1);
  }


  //////////////////////////////////////////////////////////////

  // MATRIX MORE!


  @Override
  public void resetMatrix() {
    transformMatrix.reset();
  }

  @Override
  public void applyMatrix(PMatrix2D source) {
    applyMatrixImpl(source.m00, source.m01, 0, source.m02,
                    source.m10, source.m11, 0, source.m12,
                             0,          0, 1, 0,
                             0,          0, 0, 1);
  }

  @Override
  public void applyMatrix(float n00, float n01, float n02,
                          float n10, float n11, float n12) {
    applyMatrixImpl(n00, n01, 0, n02,
                    n10, n11, 0, n12,
                      0,   0, 1,   0,
                      0,   0, 0,   1);
  }


  @Override
  public void applyMatrix(PMatrix3D source) {
    applyMatrixImpl(source.m00, source.m01, source.m02, source.m03,
                    source.m10, source.m11, source.m12, source.m13,
                    source.m20, source.m21, source.m22, source.m23,
                    source.m30, source.m31, source.m32, source.m33);
  }


  /**
   * Apply a 4x4 transformation matrix to the modelview stack.
   */
  @Override
  public void applyMatrix(float n00, float n01, float n02, float n03,
                          float n10, float n11, float n12, float n13,
                          float n20, float n21, float n22, float n23,
                          float n30, float n31, float n32, float n33) {
    applyMatrixImpl(n00, n01, n02, n03,
                    n10, n11, n12, n13,
                    n20, n21, n22, n23,
                    n30, n31, n32, n33);
  }


  protected void applyMatrixImpl(float n00, float n01, float n02, float n03,
                                 float n10, float n11, float n12, float n13,
                                 float n20, float n21, float n22, float n23,
                                 float n30, float n31, float n32, float n33) {
    transformMatrix.apply(n00, n01, n02, n03,
                    n10, n11, n12, n13,
                    n20, n21, n22, n23,
                    n30, n31, n32, n33);

  }
  //CAN'T USE PGL because Java permissions are silly
   /** Machine Epsilon for float precision. **/
  protected static float FLOAT_EPS = Float.MIN_VALUE;
  // Calculation of the Machine Epsilon for float precision. From:
  // http://en.wikipedia.org/wiki/Machine_epsilon#Approximation_using_Java
  static {
    float eps = 1.0f;

    do {
      eps /= 2.0f;
    } while ((float)(1.0 + (eps / 2.0)) != 1.0);

    FLOAT_EPS = eps;
  }

  protected static boolean same(float a, float b) {
    return Math.abs(a - b) < FLOAT_EPS;
  }


  protected static boolean diff(float a, float b) {
    return FLOAT_EPS <= Math.abs(a - b);
  }


  protected static boolean zero(float a) {
    return Math.abs(a) < FLOAT_EPS;
  }


  protected static boolean nonZero(float a) {
    return FLOAT_EPS <= Math.abs(a);
  }
  
  //end Matrix stuff from PGraphisOpenGL
  
  public boolean is3D() {
	return true;
  }
  
  public boolean is2D() {
	return false;
  }
  
  public void setColor(boolean c) {
	colorFlag = c;
	//else
	// showWarning("Cannot change color mode after beginDraw()");
  }
  
  public void setTriangleRes(int res) {
	TRIANGLE_RES = res;
	RECT_RES = res+5;
  }
  
  //override me
  protected void writeTexCoord(float u, float v) {
  
  }
  
  //override me
  protected void writeHeader() {
  }

  //override me
  protected void writeFooter() {
  }
  
  //override me
  protected void writeVertices() {
  }
  
  //override me
  protected void writeLines() {
  }

  //override me
  protected void writeFaces() {
  }
    
}
